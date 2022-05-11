# computebeta_PPT
Computes the density matrix and the optimal value for the SDP involving the PPT relaxation

## Syntax
``[RHO_PPT, BETA_PPT] = computebeta_PPT(PI)``

## Argument description
#### Input arguments
- ``PI``: input matrix

#### Output arguments
- ``RHO_PPT``: density matrix of the PPT SDP
- ``BETA_PPT``: optimal value of the PPT SDP

## Example
    >> [rho_PPT, beta_PPT] = computebeta_PPT(genMat(4, "Herm", 0))

    rho_PPT =

        0.5064 + 0.0000i  -0.3363 + 0.1808i  -0.2343 + 0.1074i   0.1173 - 0.1549i
        -0.3363 - 0.1808i   0.2878 + 0.0000i   0.1939 + 0.0123i  -0.1332 + 0.0610i
        -0.2343 - 0.1074i   0.1939 - 0.0123i   0.1312 + 0.0000i  -0.0871 + 0.0468i
        0.1173 + 0.1549i  -0.1332 - 0.0610i  -0.0871 - 0.0468i   0.0746 + 0.0000i

    beta_PPT =

        0.8318

    >> [rho_PPT, beta_PPT] = computebeta_PPT(genMat(4, "Pos", 0))

    rho_PPT =

        0.1136 + 0.0000i  -0.2141 + 0.1074i  -0.0863 - 0.0227i   0.1840 - 0.0388i
        -0.2141 - 0.1074i   0.5051 + 0.0000i   0.1411 + 0.1243i  -0.3834 - 0.1009i
        -0.0863 + 0.0227i   0.1411 - 0.1243i   0.0700 + 0.0000i  -0.1319 + 0.0662i
        0.1840 + 0.0388i  -0.3834 + 0.1009i  -0.1319 - 0.0662i   0.3112 + 0.0000i

    beta_PPT =

        0.7051

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/UpperBounds/computebeta_PPT.m)

## References
1. Peres, A. (1996). Separability criterion for density matrices. Physical Review Letters, 77(8), 1413.