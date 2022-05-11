# computebeta_prime_k
Computes the density matrix and the optimal value from the SDP involving the Bosonic Extensions relaxation

## Syntax
``[RHO_PRIME_K, BETA_PRIME_K] = computebeta_prime_k(PI, DIM_B, K, PPT)``

## Argument description
#### Input arguments
- ``PI``: input matrix
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``K``: integer that describes the number of Bosonic extensions
- ``PPT``: 1 to impose PPT criterion
           0 otherwise

#### Output arguments
- ``RHO_PRIME_K``: density matrix of the Bosonic Extensions SDP
- ``BETA_PRIME_K``: optimal value of the Bosonic Extensions SDP

## Example
    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k(genMat(2, "Herm", 0), 2, 2, 1)

    rho_prime_k =

        0.1409 + 0.0000i   0.2337 + 0.2578i
        0.2337 - 0.2578i   0.8591 + 0.0000i

    beta_prime_k =

        0.8169

    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k(genMat(2, "Pos", 0), 2, 2, 1)

    rho_prime_k =

        0.8771 + 0.0000i   0.3000 + 0.1334i
        0.3000 - 0.1334i   0.1229 + 0.0000i

    beta_prime_k =

        0.4919

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/UpperBounds/computebeta_prime_k.m)

## References
1. Chen, J., Ji, Z., Kribs, D., Lütkenhaus, N., & Zeng, B. (2014). Symmetric extension of two-qubit states. Physical Review A, 90(3), 032318.