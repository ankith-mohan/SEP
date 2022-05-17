# computegamma_rand
Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points

## Syntax
``[SIGMA_A_RAND_BEST, SIGMA_B_RAND_BEST, GAMMA_RAND_BEST, I_BEST, GAMMA_RAND_LIST] = computegamma_rand(K_LIST, L_LIST, N_SUMMANDS, DIM_A, DIM_B, THRESHOLD, N_SEESAW, N_RAND, SS_TYPE, PLUS_TWO)``

## Argument description
#### Input arguments
- ``K_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_A x DIM_A``
- ``L_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_B x DIM_B``
- ``N_SUMMANDS``: number of matrices that have to be summed over
- ``DIM_A``: dimension of Alice's subsystem
- ``DIM_B``: dimension of Bob's subsystem
- ``THRESHOLD``: threshold to determine convergence
- ``N_SEESAW``: maximum number of see-saw steps
- ``N_RAND``: number of random starting points
- ``SS_TYPE``: 
    - ``"simple"`` if ``COMPUTEGAMMA_PART``
    - ``"rev_1"`` if ``COMPUTEGAMMA_PART_REV_1``
    - ``"rev_2"`` if ``COMPUTEGAMMA_PART_REV_2``
- ``PLUS_TWO``: 
    - 2 if both the maximally mixed state and the state of uniform superposition must be added to the list of starting points
    - 0 otherwise

#### Output arguments
- ``SIGMA_A_RAND_BEST``: density matrix for Alice's subsystem
- ``SIGMA_B_RAND_BEST``: density matrix for Bob's subsystem
- ``GAMMA_RAND_BEST``: best optimal value of all the starting points
- ``I_BEST``: number of iterations that the see-saw algorithm corresponding to the best optimal value ran for
- ``GAMMA_RAND_LIST``: optimal values for all the starting points

## Example
    >> genOtimesHouseholder("test_Herm.mat", 2, 2, 2, "Herm")
    >> load("test_Herm.mat", "K_list", "L_list")
    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand_sum(K_list, L_list, 2, 2, 2, 1e-13, 100, 100, "simple", 2)

    sigma_A_rand_best =

        0.5657 + 0.0000i   0.3309 - 0.3690i
        0.3309 + 0.3690i   0.4343 + 0.0000i

    sigma_B_rand_best =

        0.0113 + 0.0000i  -0.0148 + 0.1049i
       -0.0148 - 0.1049i   0.9887 + 0.0000i

    gamma_rand_best =

        0.8669

    i_best =

        7

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand_sum(K_list, L_list, 2, 2, 2, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand_best =

        0.5657 + 0.0000i   0.3309 - 0.3690i
        0.3309 + 0.3690i   0.4343 + 0.0000i

    sigma_B_rand_best =

        0.0113 + 0.0000i  -0.0148 + 0.1049i
       -0.0148 - 0.1049i   0.9887 + 0.0000i

    gamma_rand_best =

        0.8669

    i_best =

        7

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand_sum(K_list, L_list, 2, 2, 2, 1e-13, 100, 100, "rev_2", 2)

    sigma_A_rand_best =

        0.5657 + 0.0000i   0.3309 - 0.3690i
        0.3309 + 0.3690i   0.4343 + 0.0000i

    sigma_B_rand_best =

        0.0113 + 0.0000i  -0.0148 + 0.1049i
       -0.0148 - 0.1049i   0.9887 + 0.0000i

    gamma_rand_best =

        0.8669

    i_best =

        7

    >> genOtimesHouseholder("test_Pos.mat", 2, 2, 2, "Pos")
    >> load("test_Pos.mat", "K_list", "L_list")
    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand_sum(K_list, L_list, 2, 2, 2, 1e-13, 100, 100, "simple", 2)

    sigma_A_rand_best =

        0.4422 + 0.0000i   0.4684 + 0.1651i
        0.4684 - 0.1651i   0.5578 + 0.0000i

    sigma_B_rand_best =

        0.1022 + 0.0000i  -0.2368 + 0.1889i
       -0.2368 - 0.1889i   0.8978 + 0.0000i

    gamma_rand_best =

        0.7142

    i_best =

        4

    >> genOtimesHouseholder("test_Pos.mat", 2, 2, 2, "Pos")
    >> load("test_Pos.mat", "K_list", "L_list")
    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand_sum(K_list, L_list, 2, 2, 2, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand_best =

        0.4422 + 0.0000i   0.4684 + 0.1651i
        0.4684 - 0.1651i   0.5578 + 0.0000i

    sigma_B_rand_best =

        0.1022 + 0.0000i  -0.2368 + 0.1889i
       -0.2368 - 0.1889i   0.8978 + 0.0000i

    gamma_rand_best =

        0.7142

    i_best =

        4

    >> genOtimesHouseholder("test_Pos.mat", 2, 2, 2, "Pos")
    >> load("test_Pos.mat", "K_list", "L_list")
    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand_sum(K_list, L_list, 2, 2, 2, 1e-13, 100, 100, "rev_2", 2)

    sigma_A_rand_best =

        0.4422 + 0.0000i   0.4684 + 0.1651i
        0.4684 - 0.1651i   0.5578 + 0.0000i

    sigma_B_rand_best =

        0.1022 + 0.0000i  -0.2368 + 0.1889i
       -0.2368 - 0.1889i   0.8978 + 0.0000i

    gamma_rand_best =

        0.7142

    i_best =

        4

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/LowerBounds/sum/computegamma_rand_sum.m)