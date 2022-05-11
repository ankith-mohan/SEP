# computebeta_prime_kr
Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxation

## Syntax
``[RHO_PRIME_KR, BETA_PRIME_KR] = computebeta_prime_kr(PI, DIM_B, K, PPT)``

## Argument description
#### Input arguments
- ``PI``: input matrix
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``K``: integer that describes the number of bosonic extensions
- ``PPT``: 1 to impose PPT criterion
           0 otherwise

#### Output arguments
- ``RHO_PRIME_KR``: density matrix of the Bosonic Extensions with Realignment SDP
- ``BETA_PRIME_KR``: optimal value of the Bosonic Extensions with Realignment SDP

## Example
    >> [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr(genMat(4, "Herm", 0), 4, 2, 1)

    rho_prime_kr =

        0.1454 + 0.0000i  -0.2401 - 0.0672i  -0.1256 + 0.0006i   0.2076 + 0.0570i
        -0.2401 + 0.0672i   0.4274 + 0.0000i   0.2070 - 0.0591i  -0.3691 + 0.0019i
        -0.1256 - 0.0006i   0.2070 + 0.0591i   0.1084 + 0.0000i  -0.1790 - 0.0501i
        0.2076 - 0.0570i  -0.3691 - 0.0019i  -0.1790 + 0.0501i   0.3187 + 0.0000i

    beta_prime_kr =

        0.8833

    >> [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr(genMat(4, "Pos", 0), 4, 2, 1)

    rho_prime_kr =

        0.3443 + 0.0000i  -0.4577 - 0.0095i  -0.0749 + 0.0158i   0.1000 - 0.0189i
        -0.4577 + 0.0095i   0.6086 + 0.0000i   0.0991 - 0.0230i  -0.1324 + 0.0279i
        -0.0749 - 0.0158i   0.0991 + 0.0230i   0.0170 + 0.0000i  -0.0226 - 0.0005i
        0.1000 + 0.0189i  -0.1324 - 0.0279i  -0.0226 + 0.0005i   0.0301 + 0.0000i

    beta_prime_kr =

        0.9978

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/UpperBounds/computebeta_prime_kr.m)

## References
1. Chen, J., Ji, Z., Kribs, D., Lütkenhaus, N., & Zeng, B. (2014). Symmetric extension of two-qubit states. Physical Review A, 90(3), 032318.
2. Chen, K., & Wu, L. A. (2002). A matrix realignment method for recognizing entanglement. arXiv preprint quant-ph/0205017.