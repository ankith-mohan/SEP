# computebeta_kr
Computes the density matrix and the optimal value for the SDP involving both the Symmetric Extensions as well as the Realignment relaxations

## Syntax
``[RHO_KR, BETA_KR] = computebeta_kr(PI, DIM_B, K, PPT)``

# Argument description
#### Input arguments
- ``PI``: input matrix
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``K``: integer that describes the number of symmetric extensions
- ``PPT``: 1 to impose PPT criterion
           0 otherwise

#### Output arguments
- ``RHO_KR``: density matrix of the Symmetric Extensions with Realignment SDP
- ``BETA_KR``: optimal values of the Symmetric Extensions with Realignment SDP

## Example
    >> [rho_kr, beta_kr] = computebeta_kr(genMat(4, "Herm", 0), 4, 2, 1)

    rho_kr =

        0.4066 + 0.0000i  -0.2202 + 0.1837i  -0.3146 + 0.0853i   0.1318 - 0.1883i
        -0.2202 - 0.1837i   0.2022 + 0.0000i   0.2089 + 0.0959i  -0.1564 + 0.0424i
        -0.3146 - 0.0853i   0.2089 - 0.0959i   0.2612 + 0.0000i  -0.1414 + 0.1180i
        0.1318 + 0.1883i  -0.1564 - 0.0424i  -0.1414 - 0.1180i   0.1299 + 0.0000i

    beta_kr =

        0.6325

    >> [rho_kr, beta_kr] = computebeta_kr(genMat(4, "Pos", 0), 4, 2, 1)

    rho_kr =

        0.0227 + 0.0000i   0.0250 - 0.0122i  -0.0922 - 0.0105i  -0.1068 + 0.0378i
        0.0250 + 0.0122i   0.0339 + 0.0000i  -0.0956 - 0.0608i  -0.1374 - 0.0156i
        -0.0922 + 0.0105i  -0.0956 + 0.0608i   0.3789 + 0.0000i   0.4158 - 0.2025i
        -0.1068 - 0.0378i  -0.1374 + 0.0156i   0.4158 + 0.2025i   0.5645 + 0.0000i

    beta_kr =

        0.9003

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/helpers/computebeta_kr.m)

## References
1. Chen, J., Ji, Z., Kribs, D., Lütkenhaus, N., & Zeng, B. (2014). Symmetric extension of two-qubit states. Physical Review A, 90(3), 032318.
2. Chen, K., & Wu, L. A. (2002). A matrix realignment method for recognizing entanglement. arXiv preprint quant-ph/0205017.