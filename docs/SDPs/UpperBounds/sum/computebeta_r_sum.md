# computebeta_r_sum
Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation

## Syntax
``[RHO_R, BETA_R] = computebeta_r_sum(K_LIST, L_LIST, N_SUMMANDS, DIM_B, K_PPT)``

## Argument description
#### Input arguments
- ``K_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_A x DIM_A``
- ``L_LIST``: cell of length ``N_SUMMANDS`` where each entry is a matrix of size ``DIM_B x DIM_B``
- ``N_SUMMANDS``: number of matrices that have to be summed over
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``K``: integer that describes the number of symmetric extensions
- ``PPT``: 1 to impose PPT criterion
           0 otherwise

#### Output arguments
- ``RHO_R``: density matrix of the Realignment SDP
- ``BETA_R``: optimal value of the Realignment SDP

## Example
    >> genOtimesHouseholder("test_Herm.mat", 2, 2, 2, "Herm")
    >> load("test_Herm.mat", "K_list", "L_list")
    >> [rho_r, beta_r] = computebeta_r_sum(K_list, L_list, 2, 2, 2, 1)

    rho_r =

    beta_r =

    >> genOtimesHouseholder("test_Pos.mat", 2, 2, 2, "Herm")
    >> load("test_Herm.mat", "K_list", "L_list")
    >> [rho_r, beta_r] = computebeta_r_sum(K_list, L_list, 2, 2, 2, 1)

    rho_r =

    beta_r =

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/SDPs/UpperBounds/sum/computebeta_r_sum.m)

## References
1. Chen, K., & Wu, L. A. (2002). A matrix realignment method for recognizing entanglement. arXiv preprint quant-ph/0205017.