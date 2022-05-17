#nt# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 10:51:26 2020


@author: Tobias Haug, tobiasxhaug@gmail.com
Program for NISQ SDP
For finding ground state of Hamiltonians as well as finding optimal POVMs for 
state discrimination
Support various hamiltonians and initial states.
"""
################################################################################
# Libraries
################################################################################
import os
import time
from functools import partial, reduce
import pickle

import numpy as np
import scipy
import scipy.optimize as so
from scipy.io import savemat

import qutip as qt

from helper_tools import *
################################################################################
# Options
################################################################################
optionsQt=qt.Options()
optionsQt.store_final_state=True
optionsQt.store_states=True
optionsQt.nsteps=10000

os.environ["MKL_NUM_THREADS"] = "1"
################################################################################
def genFockOp(op, position, size, levels, opdim = 0):
    opList = [qt.qeye(levels) for x in range(size - opdim)]
    opList[position] = op
    return qt.tensor(opList)


def multiply_paulis(curr_paulis, to_mult_paulis, curr_val_expansion = [],
                    to_val_expansion = []):
    """
        multiplies paulis together
    """
    new_paulis = []
    new_val = []
    for i in range(len(curr_paulis)):
        for j in range(len(to_mult_paulis)):
            add_pauli = np.zeros(len(curr_paulis[i]), dtype = int)

            for k in range(len(curr_paulis[i])):
                if (curr_paulis[i][k] == 1 and to_mult_paulis[j][k] == 2) or \
                    (curr_paulis[i][k] == 2 and to_mult_paulis[j][k] == 1):
                    add_pauli[k] = 3
                else:
                    add_pauli[k] = np.abs(curr_paulis[i][k] \
                                            - to_mult_paulis[j][k])

            new_paulis.append(add_pauli)
            if len(curr_val_expansion) > 0:
                new_val.append(curr_val_expansion[i] * to_val_expansion[j])

    # new_paulis = list(np.unique(new_paulis, axis = 0))
    new_paulis, inverse_array = list(np.unique(new_paulis, axis = 0, \
                                    return_inverse = True))
    new_val = np.array(new_val)
    new_val_unique = np.zeros(len(new_paulis))
    
    # reconstruct weight values for each pauli, 
    # for paulis which occur multiple times are added up
    for i in range(len(new_val)):
        new_val_unique[inverse_array[i]] += np.abs(new_val[i])
    
    return new_paulis, new_val_unique


def get_ini_state(ini_state_type, levels, n_qubits, depth, opcsign, annealtime,
                    model, Hvalues, J, H):
    """
        Initialize state
        Args:
            - ini_type
                state to use as ansatz for k-moment expansion 
                initial state for ansatz state via K moment
                Options:
                    - 0: + product state
                    - 1: zero state
                    - 2: random circuit
                    - 15: circuit annealing with fixed anneal time
                    - 16: circuit annealing where annealing time is optimized
            - model
                model Hamiltonian used
                Options:
                    - 0: transverse ising
                    - 1: Heisenberg
                    - 14: transverse ising model with longitudinal field
                    - 30: POVM
    """
    global anneal_time_opt
    #get initial state

    if ini_state_type == 0: # product state plus state
        initial_state = qt.tensor([qt.basis(levels, 1) \
                                + qt.basis(levels, 0) for _ in range(n_qubits)])
        # initial_state=qt.tensor([qt.basis(levels, 1)] \
        #                       + [qt.basis(levels, 0) for _ in range(L-1)]) 
        # tjis was used for paper to compare against imag time evolution
    
    elif ini_state_type == 1: # all 0
        initial_state = qt.tensor([qt.basis(levels, 0) for _ in range(n_qubits)])
    
    elif ini_state_type == 2: # random state
        rand_angles = np.random.rand(depth, n_qubits) * 2 * np.pi
        rand_pauli = np.random.randint(1, 4, [depth, n_qubits])

        entangling_layer = prod([opcsign[j] for j in range(n_qubits - 1)])
        initial_state = qt.tensor([qt.basis(levels, 0) for _ in range(n_qubits)])
        initial_state = qt.tensor([qt.qip.operations.ry(np.pi/4) \
                                    for _ in range(n_qubits)]) * initial_state
        
        for j in range(depth):
            rot_op = []
            for k in range(n_qubits):
                angle = rand_angles[j][k]
                if rand_pauli[j][k] == 1:
                    rot_op.append(qt.qip.operations.rx(angle))
                elif rand_pauli[j][k] == 2:
                    rot_op.append(qt.qip.operations.ry(angle))
                elif rand_pauli[j][k] == 3:
                    rot_op.append(qt.qip.operations.rz(angle))

            initial_state = qt.tensor(rot_op) * initial_state
            initial_state = entangling_layer * initial_state

    elif ini_state_type == 15: # time annealing ansatz
        initial_state = qt.tensor([qt.basis(levels, 1) - qt.basis(levels, 0) \
                                    for _ in range(n_qubits)])
        initial_state /= initial_state.norm()
        if annealtime != 0:
            ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
        if model == 0:
            # H1 = 0 # ZZ terms
            # for i in range(len(Hstrings)//2):
            #     H1 += Hvalues[i] * get_pauli_op(Hstrings[i])
            # H2 = 0 # X terms
            # for i in range(len(Hstrings)//2):
            #     H2 += Hvalues[i + len(Hstrings)//2] \
            #           * get_pauli_op(Hstrings[i + len(Hstrings)//2])
                
            # initial_state = qt.tensor([qt.basis(levels, 1) \
            #                             - qt.basis(levels, 0) \
            #                             for _ in range(n_qubits)])
            # initial_state /= initial_state.norm()
            if annealtime != 0:
                # ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
                for i in range(depth):
                    ZZrot = (-1j * Hvalues[0] * annealtime / depth * ZZ \
                            * (i+1)).expm()
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_2toN(ZZrot, \
                                            n_qubits, j, (j+1) % n_qubits) \
                                            * initial_state
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rx(2 \
                                            * Hvalues[-1] * annealtime), \
                                            n_qubits, j) * initial_state
        elif model == 14: ## transverse ising with longitudinal field
            # initial_state = qt.tensor([qt.basis(levels, 1) \
            #                             - qt.basis(levels, 0) \
            #                             for _ in range(n_qubits)])
            # initial_state /= initial_state.norm()
            if annealtime != 0:
                # ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
                for i in range(depth):
                    for j in range(n_qubits):
                        ZZrot = (-1j * J * Hvalues[j] * annealtime / depth \
                                * ZZ * (i+1)).expm()
                        initial_state = qt.qip.operations.gate_expand_2toN(ZZrot, \
                                            n_qubits, j, (j+1) % n_qubits) \
                                            * initial_state
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rz(2 \
                                            * Hvalues[j + 2 * n_qubits] \
                                            * annealtime / depth * (i+1)), \
                                            n_qubits, j) * initial_state
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rx(2 \
                                            * Hvalues[j + n_qubits] \
                                            * annealtime), n_qubits, j) \
                                            * initial_state
                     
    elif ini_state_type == 16: ## optimize annealing time
        def get_energy_ramp(x):
            annealtime = x[0]
            initial_state = qt.tensor([qt.basis(levels, 1) \
                                        - qt.basis(levels, 0) \
                                        for _ in range(n_qubits)])
            initial_state /= initial_state.norm()
            if annealtime != 0:
                ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
            if model == 0:
                if annealtime != 0:
                    # ZZ=qt.tensor([qt.sigmaz(),qt.sigmaz()])
                    for i in range(depth):
                        ZZrot = (-1j * Hvalues[0] * annealtime / depth * ZZ \
                                    * (i+1)).expm()
                        for j in range(n_qubits):
                            initial_state = qt.qip.operations.gate_expand_2toN(ZZrot, \
                                                n_qubits, j, (j+1) % n_qubits) \
                                                * initial_state
                        for j in range(n_qubits):
                            initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rx(2 \
                                                * Hvalues[-1] * annealtime), \
                                                n_qubits, j) * initial_state            
            elif model == 14: ## transverse ising with longitudinal field
                if annealtime != 0:
                    # ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
                    for i in range(depth):
                        for j in range(n_qubits):
                            ZZrot = (-1j * J * Hvalues[j] * annealtime / depth \
                                    * ZZ * (i+1)).expm()
                            initial_state = qt.qip.operations.gate_expand_2toN(ZZrot, \
                                                n_qubits, j, (j+1) % n_qubits) \
                                                * initial_state
                        for j in range(n_qubits):
                            initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rz(2 \
                                                * Hvalues[j + 2 * n_qubits] \
                                                * annealtime / depth * (i+1)), \
                                                n_qubits, j) * initial_state
                        for j in range(n_qubits):
                            initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rx(2 \
                                                * Hvalues[j + n_qubits] \
                                                * annealtime), n_qubits, j) \
                                                * initial_state
            circuit_energy = qt.expect(H, initial_state)
            return circuit_energy
        
        annealtime_ini = annealtime
        ramp_method = "Nelder-Mead"
        options = {"maxiter": 20}
        res = so.minimize(get_energy_ramp, [annealtime_ini], \
                            method = ramp_method, options = options)
        anneal_time_opt = res["x"][0]
        print(f"anneal time found: {anneal_time_opt}", res["fun"], res["nit"])

        initial_state = qt.tensor([qt.basis(levels, 1) - qt.basis(levels, 0) \
                                    for _ in range(n_qubits)])
        initial_state /= initial_state.norm()
        annealtime_run = anneal_time_opt
        if model == 0:
            if annealtime_run != 0:
                ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
                for i in range(depth):
                    ZZrot = (-1j * Hvalues[0] * annealtime_run / depth * ZZ \
                                           * (i+1)).expm()
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_2toN(ZZrot, \
                                            n_qubits, j, (j+1) % n_qubits) \
                                            * initial_state
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rx(2 \
                                            * Hvalues[-1] * annealtime_run), \
                                            n_qubits, j) * initial_state
            
        elif model == 14: ## transverse ising  with longitudinal field
            if annealtime_run != 0:
                ZZ = qt.tensor([qt.sigmaz(), qt.sigmaz()])
                for i in range(depth):
                    for j in range(n_qubits):
                        ZZrot = (-1j * J * Hvalues[j] * annealtime_run / depth \
                                    * ZZ * (i+1)).expm()
                        initial_state = qt.qip.operations.gate_expand_2toN(ZZrot, \
                                            n_qubits, j, (j+1) % n_qubits) \
                                            * initial_state
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rz(2 \
                                            * Hvalues[j + 2 * n_qubits] \
                                            * annealtime_run / depth * (i+1)), \
                                            n_qubits, j) * initial_state
                    for j in range(n_qubits):
                        initial_state = qt.qip.operations.gate_expand_1toN(qt.qip.operations.rx(2 \
                                            * Hvalues[j + n_qubits] \
                                            * annealtime_run), n_qubits, j) \
                                            * initial_state
    
    initial_state/=initial_state.norm()
    return initial_state


def get_Hamiltonian_string(L, model, J, h, g):
    """
        get pauli strings for models
        Args:
            - J: J parameter
            - h: h parameter
    """
    Hstrings, Hvalues = [], []
    Hoffset = 0
    
    if model == 0: # ising
        if J != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 3
                paulistring[(i+1) % L] = 3
                Hstrings.append(list(paulistring))
                Hvalues.append(0.5 * J)
        if h != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 1
                Hstrings.append(list(paulistring))
                Hvalues.append(0.5 * h)
    
    elif model == 1: # heisenberg
        if h != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 3
                paulistring[(i+1) % L] = 3
                Hstrings.append(list(paulistring))
                Hvalues.append(h)
        if J != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 2
                paulistring[(i+1) % L] = 2
                Hstrings.append(list(paulistring))
                Hvalues.append(J)
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 1
                paulistring[(i+1) % L] = 1
                Hstrings.append(list(paulistring))
                Hvalues.append(J)
            
    elif model == 14: # transverse ising with longitudinal field
        if J != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 3
                paulistring[(i+1) % L] = 3
                Hstrings.append(list(paulistring))
                Hvalues.append(J)
        if h != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 1
                Hstrings.append(list(paulistring))
                Hvalues.append(h)     
        if g != 0:
            for i in range(L):
                paulistring = np.zeros(L, dtype = int)
                paulistring[i] = 3
                Hstrings.append(list(paulistring))
                Hvalues.append(g)
    
    HpauliFactor = np.zeros([len(Hstrings), L, 4])
    for i in range(len(Hstrings)):
        pauliFactor = np.zeros([L, 4])
        for j in range(L):
            if Hstrings[i][j] == 0:
                pauliFactor[j] = [1,0,0,0]
            elif Hstrings[i][j] == 1:
                pauliFactor[j] = [0,1,0,0]
            elif Hstrings[i][j] == 2:
                pauliFactor[j] = [0,0,1,0]
            elif Hstrings[i][j] == 3:
                pauliFactor[j] = [0,0,0,1]
        HpauliFactor[i] = pauliFactor
    
    return Hstrings, HpauliFactor, Hvalues, Hoffset


def get_pauli_op(pauli_string, opId, opX, opY, opZ):
    """
        make operator from pauli string
    """
    pauli_circuit = opId
    for i in range(len(pauli_string)):
        if pauli_string[i] != 0:
            if pauli_string[i] == 1:
                pauli_circuit = pauli_circuit * opX[i]
            elif pauli_string[i] == 2:
                pauli_circuit = pauli_circuit * opY[i]
            elif pauli_string[i] == 3:
                pauli_circuit = pauli_circuit * opZ[i]
    return pauli_circuit


def decomposePauli(H):
    """
        Decompose Hermitian matrix H into Pauli matrices
    """
    S = get_pauli_mats()
    dim_matrix = np.shape(H)[0]
    n_qubits = int(np.log2(dim_matrix))
    assert dim_matrix == 2**n_qubits, IOError("Matrix is not power of 2!")
    hilbertspace = 2**n_qubits
    n_paulis = 4**n_qubits
    pauli_list = np.zeros([n_paulis, n_qubits], dtype = int)
    for k in range(n_paulis):
        pauli_list[k, :] = numberToBase(k, 4, n_qubits)
    weights = np.zeros(n_paulis, dtype = np.complex128)
    for k in range(n_paulis):
        pauli = S[pauli_list[k][0]]
        for n in range(1, n_qubits):
            pauli = np.kron(pauli, S[pauli_list[k][n]])
        # weights[k] = 1 / hilbertspace * (H.conjugate().transpose() @ pauli).trace()
        weights[k] = 1 / hilbertspace * (pauli @ H).trace()
    return pauli_list, weights


def reconstructMatrix(pauli_list, weights):
    S = get_pauli_mats()
    n_paulis = len(weights)
    n_qubits = np.shape(pauli_list)[1]
    H = np.zeros([2**n_qubits, 2**n_qubits], dtype = np.complex128)
    for k in range(n_paulis):
        if weights[k] != 0:
            pauli = S[pauli_list[k][0]]
            for n in range(1,n_qubits):
                pauli = np.kron(pauli, S[pauli_list[k][n]])
            H += weights[k] * pauli
    return H