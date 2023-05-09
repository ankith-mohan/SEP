# List of functions

## Upper bounds
These are the functions that correspond to relaxations of the separability set.
- [computebeta_PPT](/SDPs/UpperBounds/operator/computebeta_PPT.md): Computes the density matrix and the optimal value for the SDP involving PPT relaxations.
- [computebeta_k](/SDPs/UpperBounds/operator/computebeta_k.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation.
- [computebeta_kr](/SDPs/UpperBounds/operator/computebeta_kr.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations.
- [computebeta_prime_k](/SDPs/UpperBounds/operator/computebeta_prime_k.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation.
- [computebeta_prime_kr](/SDPs/UpperBounds/operator/computebeta_prime_kr.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxations.
- [computebeta_r](/SDPs/UpperBounds/operator/computebeta_r.md): Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation.

If getDE is used to compute the overlap measurement matrices ($D$ and $E$), then the following functions must be used.
- [computebeta_PPT_DE](/SDPs/UpperBounds/operator_NISQ/computebeta_PPT_DE.md): Computes the density matrix and the optimal value for the SDP involving PPT relaxations.
- [computebeta_k_DE](/SDPs/UpperBounds/operator_NISQ/computebeta_k_DE.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation.
- [computebeta_kr_DE](/SDPs/UpperBounds/operator_NISQ/computebeta_kr_DE.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations.
- [computebeta_prime_k_DE](/SDPs/UpperBounds/operator_NISQ/computebeta_prime_k_DE.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation.
- [computebeta_prime_kr_DE](/SDPs/UpperBounds/operator_NISQ/computebeta_prime_kr_DE.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxations.
- [computebeta_r_DE](/SDPs/UpperBounds/operator_NISQ/computebeta_r_DE.md): Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation.

If (1) [genOtimesHouseholder](/helpers/genOtimesHouseholder.md) is used to generate the cells of matrices or (2) the matrices are special case instances, then the following counterparts must be used.
- [computebeta_PPT_sum](/SDPs/UpperBounds/sum_operator/computebeta_PPT_sum.md): Computes the density matrix and the optimal value for the SDP involving PPT relaxations.
- [computebeta_k_sum](/SDPs/UpperBounds/sum_operator/computebeta_k_sum.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation.
- [computebeta_kr_sum](/SDPs/UpperBounds/sum_operator/computebeta_kr_sum.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations.
- [computebeta_prime_k_sum](/SDPs/UpperBounds_operator/sum/computebeta_prime_k_sum.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation.
- [computebeta_prime_kr_sum](/SDPs/UpperBounds/sum_operator/computebeta_prime_kr_sum.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxations.
- [computebeta_r_sum](/SDPs/UpperBounds/sum_operator/computebeta_r_sum.md): Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation.

If (1) [genOtimesHouseholder](/helpers/genOtimesHouseholder.md) is used to generate the cells of matrices or (2) the matrices are special case instances, and getDE is used to compute the overlap measurement matrices ($D$ and $E$), then the following must be used.
- [computebeta_PPT_DE_sum](/SDPs/UpperBounds/sum_operator_NISQ/computebeta_PPT_DE_sum.md): Computes the density matrix and the optimal value for the SDP involving PPT relaxations.
- [computebeta_k_DE_sum](/SDPs/UpperBounds/sum_operator_NISQ/computebeta_k_DE_sum.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation.
- [computebeta_kr_DE_sum](/SDPs/UpperBounds/sum_operator_NISQ/computebeta_kr_DE_sum.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations.
- [computebeta_prime_k_DE_sum](/SDPs/UpperBounds/sum_operator_NISQ/computebeta_prime_k_DE_sum.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation.
- [computebeta_prime_kr_DE_sum](/SDPs/UpperBounds/sum_operator_NISQ/computebeta_prime_kr_DE_sum.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxations.
- [computebeta_r_DE_sum](/SDPs/UpperBounds/sum_operator_NISQ/computebeta_r_DE_sum.md): Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation.

## Lower bounds
These functions compute the lower bound using the See-saw method.
- [computegamma](/SDPs/LowerBounds/operator/computegamma.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for the see-saw algorithm.
- [computegamma_rev_1](/SDPs/LowerBounds/operator/computegamma_rev_1.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 1 of the see-saw algorithm.
- [computegamma_rev_2](/SDPs/LowerBounds/operator/computegamma_rev_2.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 2 of the see-saw algorithm.
- [computegamma_rand](/SDPs/LowerBounds/operator/computegamma_rand.md): Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points.

If getDE is used to compute the overlap measurement matrices ($D$ and $E$), then the following functions must be used.
- [computegamma_DE](/SDPs/LowerBounds/operator_NISQ/computegamma_DE.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for the see-saw algorithm.
- [computegamma_rev_1_DE](/SDPs/LowerBounds/operator_NISQ/computegamma_rev_1_DE.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 1 of the see-saw algorithm.
- [computegamma_rev_2_DE](/SDPs/LowerBounds/operator_NISQ/computegamma_rev_2_DE.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 2 of the see-saw algorithm.
- [computegamma_rand_DE](/SDPs/LowerBounds/operator_NISQ/computegamma_rand_DE.md): Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points.

If (1) [genOtimesHouseholder](/helpers/genOtimesHouseholder.md) is used to generate the cells of matrices or (2) the matrices are special case instances, then the following counterparts must be used.
- [computegamma_sum](/SDPs/LowerBounds/sum_operator/computegamma_sum.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for the see-saw algorithm.
- [computegamma_rev_1_sum](/SDPs/LowerBounds/sum_operator/computegamma_rev_1_sum.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 1 of the see-saw algorithm.
- [computegamma_rev_2_sum](/SDPs/LowerBounds/sum_operator/computegamma_rev_2_sum.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 2 of the see-saw algorithm.
- [computegamma_rand_sum](/SDPs/LowerBounds/sum_operator/computegamma_rand_sum.md): Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points.

If (1) [genOtimesHouseholder](/helpers/genOtimesHouseholder.md) is used to generate the cells of matrices or (2) the matrices are special case instances, and getDE is used to compute the overlap measurement matrices ($D$ and $E$), then the following must be used.
- [computegamma_DE_sum](/SDPs/LowerBounds/sum_operator_NISQ/computegamma_DE_sum.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for the see-saw algorithm.
- [computegamma_rev_1_DE_sum](/SDPs/LowerBounds/sum_operator_NISQ/computegamma_rev_1_DE_sum.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 1 of the see-saw algorithm.
- [computegamma_rev_2_DE_sum](/SDPs/LowerBounds/sum_operator_NISQ/computegamma_rev_2_DE_sum.md): Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for revision 2 of the see-saw algorithm.
- [computegamma_rand_DE_sum](/SDPs/LowerBounds/sum_operator_NISQ/computegamma_rand_DE_sum.md): Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points.


## NISQ
- [consAnsatzStates](/NISQ/consAnsatzStates.md): Constructs the ansatz states from a reference state using Krylov expansion.

## Driver functions
These are the functions that can be used to monitor computations of the bounds
- [setupExpt](/main/setupExpt.md): Computes the upper and lower bounds for ``N_OPS`` instances.
- [executeExpt](/main/executeExpt.md): Creates the necessary folder structure for ``SETUPEXPT``.
- [TIHam](/main/TIHam.md): 
- [diagHam](/main/diagHam.md):
- [diagKLHam](/main/diagKLHam.md):


## Helper functions
These are functions that serve auxillary roles with SEP. As an end-user they might not be useful, but are listed anyway.
- [HSIP](/helpers/HSIP.md): Computes the Hilbert-Schmidt inner product between the input matrices.
- [genHerm](/helpers/genHerm.md): Generates a complex Hermitian matrix of the input dimensionality.
- [genHouseholder](/helpers/genHouseholder.md): Generates a complex Hermitian or Positive Semidefinite matrix of the input dimensionality.
- [genMat](/helpers/genMat.md): Generates a complex Hermitian or Positive Semidefinite matrix of the input dimensionality normalized by its trace norm.
- [genOtimesHouseholder](/helpers/genOtimesHouseholder.md): Generates two cells comprising of a set of complex Hermitian or Positive Semidefinite matrices of the respective dimensionalities.
- [genPos](/helpers/genPos.md): Generates a Positive Semidefinite matrix of the input dimensionality.
- [genEvals](/helpers/genEvals.md): 
- [computeEvals](/SDPS/computeEvals.md): Compute the largest and smallest eigenvalues given the overlap measurement matrices ($D$ and $E$).
