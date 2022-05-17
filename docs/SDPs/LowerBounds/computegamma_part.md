# computegamma_part
Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for the see-saw algorithm

## Syntax
``[SIGMA_A, SIGMA_B, GAMMA, I] = computegamma_part(PI, SIGMA, DIM, THRESHOLD, N_SEESAW, OPT)``

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
    >> [sigma_A, sigma_B, gamma, i] = computegamma_part(genMat(4, "Herm", 0), genMat(2, "Herm", 0), 2, 1e-13, 100, "sigma_A")

    sigma_A =

        0.4952 + 0.0000i   0.4651 + 0.1835i
        0.4651 - 0.1835i   0.5048 + 0.0000i

    sigma_B =

        0.5038 + 0.0000i   0.3007 + 0.3995i
        0.3007 - 0.3995i   0.4962 + 0.0000i

    gamma =

        0.5664

    i =

        12

    >> [sigma_A, sigma_B, gamma, i] = computegamma_part(genMat(4, "Pos", 0), genMat(2, "Pos", 0), 2, 1e-13, 100, "sigma_B")

    sigma_A =

        0.7814 + 0.0000i  -0.1280 + 0.3930i
       -0.1280 - 0.3930i   0.2186 + 0.0000i


    sigma_B =

        0.1404 + 0.0000i  -0.1789 + 0.2978i
       -0.1789 - 0.2978i   0.8596 + 0.0000i

    gamma =

        0.7234

    i =

        8

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/LowerBounds/computegamma_part.m)