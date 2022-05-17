## computegamma_part_rev_2
Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 2 of the see-saw algorithm

## Syntax
``[SIGMA_A, SIGMA_B, GAMMA, I] = computegamma_part_rev_2(PI, SIGMA, DIM, THRESHOLD, N_SEESAW, OPT)``

## Argument description
#### Input arguments
- ``PI``: input matrix
- ``SIGMA``: starting matrix for the see-saw algorithm
- ``DIM``: dimension of the starting matrix
- ``THRESHOLD``: threshold to determine convergence
- ``N_SEESAW``: maximum number of see-saw steps
- ``OPT``: 
    - ``"sigma_A"`` if starting matrix is for Alice's subsystem
    - ``"sigma_B"`` if starting matrix is for Bob's subsystem

#### Output arguments
- ``SIGMA_A``: density matrix for Alice's subsystem
- ``SIGMA_B``: density matrix for Bob's subsystem
- ``GAMMA``: optimal value of see-saw algorithm
- ``I``: number of iterations that the see-saw algorithm ran

## Example
    >> [sigma_A, sigma_B, gamma, i] = computegamma_part_rev_2(genMat(4, "Herm", 0), genMat(2, "Herm", 0), 2, 1e-13, 100, "sigma_A")

    sigma_A =

        0.4972 + 0.0000i  -0.3812 - 0.3235i
       -0.3812 + 0.3235i   0.5028 + 0.0000i

    sigma_B =

        0.1100 + 0.0000i  -0.1208 - 0.2887i
       -0.1208 + 0.2887i   0.8900 + 0.0000i

    gamma =

       -0.1008

    i =

        16

    >> [sigma_A, sigma_B, gamma, i] = computegamma_part_rev_2(genMat(4, "Pos", 0), genMat(2, "Pos", 0), 2, 1e-13, 100, "sigma_B")

    sigma_A =

        0.0925 + 0.0000i   0.0696 - 0.2813i
        0.0696 + 0.2813i   0.9075 + 0.0000i

    sigma_B =

        0.4326 + 0.0000i  -0.4938 + 0.0400i
       -0.4938 - 0.0400i   0.5674 + 0.0000i

    gamma =

        0.7548

    i =

        15

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/LowerBounds/computegamma_part_rev_2.m)