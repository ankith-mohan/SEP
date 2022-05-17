# computegamma_part_rev_2_sum
Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 2 of the see-saw algorithm

## Syntax
``[SIGMA_A, SIGMA_B, GAMMA, I] = computegamma_part_rev_2_sum(K_LIST, L_LIST, N_SUMMANDS, SIGMA, DIM, THRESHOLD, N_SEESAW, OPT)``

## Argument description
#### Input arguments
- ``K_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_A x DIM_A``
- ``L_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_B x DIM_B``
- ``N_SUMMANDS``: number of matrices that have to be summed over
- ``SIGMA``: starting matrix for the see-saw algorithm
- ``DIM``: dimension of the starting matrix
- ``THRESHOLD``: threshold to determine convergence
- ``N_SEESAW``: maximum number of see-saw steps
- ``OPT``: 
    - ``"sigma_A"`` is starting matrix is for Alice's subsystem
    - ``"sigma_B"`` is starting matrix for Bob's subsystem

#### Output arguments
- ``SIGMA_A``: density matrix for Alice's subsystem
- ``SIGMA_B``: density matrix for Bob's subsystem
- ``GAMMA``: optimal value of see-saw algorithm
- ``I``: number of iterations that the see-saw algorithm ran

## Example
    >> genOtimesHouseholder("test_Herm.mat", 2, 2, 2, "Herm")
    >> load("test_Herm.mat", "K_list", "L_list")
    >> [sigma_A, sigma_B, gamma, i] = computegamma_part_rev_2_sum(K_list, L_list, 2, genMat(2, "Herm", 0), 2, 1e-13, 100, "sigma_A")

    sigma_A =

        0.5657 + 0.0000i   0.3309 - 0.3690i
        0.3309 + 0.3690i   0.4343 + 0.0000i

    sigma_B =

        0.0113 + 0.0000i  -0.0148 + 0.1049i
       -0.0148 - 0.1049i   0.9887 + 0.0000i

    gamma =

        0.8669

    i =

        8

    >> genOtimesHouseholder("test_Pos.mat", 2, 2, 2, "Pos")
    >> load("test_Pos.mat", "K_list", "L_list")
    >> [sigma_A, sigma_B, gamma, i] = computegamma_part_rev_2_sum(K_list, L_list, 2, genMat(2, "Pos", 0), 2, 1e-13, 100, "sigma_B")

    sigma_A =

        0.4422 + 0.0000i   0.4684 + 0.1651i
        0.4684 - 0.1651i   0.5578 + 0.0000i

    sigma_B =

        0.1022 + 0.0000i  -0.2368 + 0.1889i
       -0.2368 - 0.1889i   0.8978 + 0.0000i

    gamma =

        0.7142

    i =

        4

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/LowerBounds/sum/computegamma_part_rev_2_sum.m)