# computebeta_prime_k_sum
Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation

## Syntax
``[RHO_PRIME_K, BETA_PRIME_K] = computebeta_prime_k_sum(K_LIST, L_LIST, N_SUMMANDS, DIM_B, K, PPT)``

## Argument description
#### Input arguments
- ``K_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_A x DIM_A``
- ``L_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_B x DIM_B``
- ``N_SUMMANDS``: number of matrices that have to be summed over
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``PPT``: 1 to impose PPT criterion
           0 otherwise

#### Output arguments
- ``RHO_PRIME_K``: density matrix of the Bosonic Extensions SDP
- ``BETA_PRIME_K``: optimal value of the Bosonic Extensions SDP

## Example
    >> genOtimesHouseholder("test_Herm.mat", 2, 2, 2, "Herm")
    >> load("test_Herm.mat", "K_list", "L_list")
    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k_sum(K_list, L_list, 2, 2, 2, 1)

    rho_prime_k =

    beta_prime_k =

    >> genOtimesHouseholder("test_Pos.mat", 2, 2, 2, "Herm")
    >> load("test_Herm.mat", "K_list", "L_list")
    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k_sum(K_list, L_list, 2, 2, 2, 1)

    rho_prime_k =

    beta_prime_k =

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/UpperBounds/sum/computebeta_prime_k_sum.m)

## References
1. Chen, J., Ji, Z., Kribs, D., Lütkenhaus, N., & Zeng, B. (2014). Symmetric extension of two-qubit states. Physical Review A, 90(3), 032318.