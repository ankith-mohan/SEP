# computebeta_k
Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation

## Syntax
``[RHO_K, BETA_K] = computebeta_k(PI, DIM_B, K, PPT)``

## Argument description
#### Input arguments
- ``PI``: input matrix
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``K``: integer that describes the number of symmetric extensions
- ``PPT``: 1 to impose PPT criterion
           0 otherwise

#### Output arguments
- ``RHO_K``: density matrix of the Symmetric Extensions SDP
- ``BETA_K``: optimal value of the Symmetric Extensions SDP

## Example
    >> [rho_k, beta_k] = computebeta_k(genMat(2, "Herm", 0), 2, 2, 1)

    rho_k = 

        0.8740 + 0.0000i   0.3102 + 0.1179i
        0.3102 - 0.1179i   0.1260 + 0.0000i

    beta_k = 

        -0.5315

    >> [rho_k, beta_k] = computebeta_k(genMat(2, "Pos", 0), 2, 2, 1)

    rho_k =

        0.5844 + 0.0000i   0.2183 - 0.4419i
        0.2183 + 0.4419i   0.4156 + 0.0000i

    beta_k =

        0.1088

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/UpperBounds/computebeta_k.m)

## References
1. Chen, J., Ji, Z., Kribs, D., Lütkenhaus, N., & Zeng, B. (2014). Symmetric extension of two-qubit states. Physical Review A, 90(3), 032318.