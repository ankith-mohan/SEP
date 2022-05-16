# computegamma_rand
Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points

## Syntax
``[SIGMA_A_RAND_BEST, SIGMA_B_RAND_BEST, GAMMA_RAND_BEST, I_BEST, GAMMA_RAND_LIST] = computegamma_rand(PI, DIM_A, DIM_B, THRESHOLD, N_SEESAW, N_RAND, SS_TYPE, PLUS_TWO)``

## Argument description
#### Input arguments
- ``PI``: input matrix
- ``DIM_A``: dimension of Alice's subsystem
- ``DIM_B``: dimension of Bob's subsystem
- ``THRESHOLD``: threshold to determine convergence
- ``N_SEESAW``: maximum number of see-saw steps
- ``N_RAND``: number of random starting points
- ``SS_TYPE``: ``"simple"`` if ``COMPUTEGAMMA_PART``
               ``"rev_1"`` if ``COMPUTEGAMMA_PART_REV_1``
               ``"rev_2"`` if ``COMPUTEGAMMA_PART_REV_2``
- ``PLUS_TWO``: 2 if both the maximally mixed state and the state of uniform superposition must be added to the list of starting points
                0 otherwise

#### Output arguments
- ``SIGMA_A_RAND_BEST``: density matrix for Alice's subsystem
- ``SIGMA_B_RAND_BEST``: density matrix for Bob's subsystem
- ``GAMMA_RAND_BEST``: best optimal value of all the starting points
- ``I_BEST``: number of iterations that the see-saw algorithm corresponding to the best optimal value ran for
- ``GAMMA_RAND_LIST``: optimal values for all the starting points

## Example
    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand(genMat(4, "Herm", 0), 2, 2, 1e-13, 100, 100, "simple", 2)

    sigma_A_rand_best =

        0.8490 + 0.0000i   0.1390 + 0.3299i
        0.1390 - 0.3299i   0.1510 + 0.0000i

    sigma_B_rand_best =

        0.6621 + 0.0000i  -0.3973 + 0.2567i
        -0.3973 - 0.2567i   0.3379 + 0.0000i

    gamma_rand_best =

        0.0651

    i_best =

        9

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand(genMat(4, "Herm", 0), 2, 2, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand_best =

        0.8490 + 0.0000i   0.1390 + 0.3299i
        0.1390 - 0.3299i   0.1510 + 0.0000i

    sigma_B_rand_best =

        0.6621 + 0.0000i  -0.3973 + 0.2567i
       -0.3973 - 0.2567i   0.3379 + 0.0000i

    gamma_rand_best =

        0.0651

    i_best =

        12

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand(genMat(4, "Herm", 0), 2, 2, 1e-13, 100, 100, "rev_2", 2)

    sigma_A_rand_best =

        0.4972 + 0.0000i  -0.3812 - 0.3235i
       -0.3812 + 0.3235i   0.5028 + 0.0000i

    sigma_B_rand_best =

        0.1100 + 0.0000i  -0.1208 - 0.2887i
       -0.1208 + 0.2887i   0.8900 + 0.0000i

    gamma_rand_best =

        -0.1008

    i_best =

        100

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand(genMat(4, "Pos", 0), 2, 2, 1e-13, 100, 100, "simple", 2)

    sigma_A_rand_best =

        0.9506 + 0.0000i  -0.1751 + 0.1277i
       -0.1751 - 0.1277i   0.0494 + 0.0000i

    sigma_B_rand_best =

        0.9920 + 0.0000i   0.0724 - 0.0524i
        0.0724 + 0.0524i   0.0080 + 0.0000i

    gamma_rand_best =

        0.8475

    i_best =

        14

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand(genMat(4, "Pos", 0), 2, 2, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand_best =

        0.9506 + 0.0000i  -0.1751 + 0.1277i
       -0.1751 - 0.1277i   0.0494 + 0.0000i

    sigma_B_rand_best =

        0.9920 + 0.0000i   0.0724 - 0.0524i
        0.0724 + 0.0524i   0.0080 + 0.0000i

    gamma_rand_best =

        0.8475

    i_best =

        14

    >> [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, i_best, ~] = computegamma_rand(genMat(4, "Pos", 0), 2, 2, 1e-13, 100, 100, "rev_2", 2)

    sigma_A_rand_best =

        0.0925 + 0.0000i   0.0696 - 0.2813i
        0.0696 + 0.2813i   0.9075 + 0.0000i

    sigma_B_rand_best =

        0.4326 + 0.0000i  -0.4938 + 0.0400i
       -0.4938 - 0.0400i   0.5674 + 0.0000i

    gamma_rand_best =

        0.7548

    i_best =

        14

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/LowerBounds/computegamma_rand.m)