**SEP** is a MATLAB toolbox for computing bounds on the separability problem.

## List of functions
### Upper bounds
These are the functions that correspond to relaxations of the separability set.
- [computebeta_PPT](/SDPs/UpperBounds/computebeta_PPT.md): Computes the density matrix and the optimal value for the SDP involving PPT relaxations.
- [computebeta_k](/SDPs/UpperBounds/computebeta_k.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation.
- [computebeta_kr](/SDPs/UpperBounds/computebeta_kr.md): Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations.
- [computebeta_prime_k](/SDPs/UpperBounds/computebeta_prime_k.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation.
- [computebeta_prime_kr](/SDPs/UpperBounds/computebeta_prime_kr.md): Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxations.
- [computebeta_r](/SDPs/UpperBounds/computebeta_r): Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation.

### Helper functions
These are functions that serve auxillary roles with SEP. As an end-user they might not be useful, but are listed anyway.
- [HSIP](/helpers/HSIP.md): Computes the Hilbert-Schmidt inner product between the input matrices.
- [genHerm](/helpers/genHerm.md): Generates a complex Hermitian matrix of the input dimensionality.
- [genHouseholder](/helpers/genHouseholder.md): Generates a complex Hermitian or Positive Semidefinite matrix of the input dimensionality.
- [genMat](/helpers/genMat.md): Generates a complex Hermitian or Positive Semidefinite matrix of the input dimensionality normalized by its trace norm.
- [genOtimesHouseholder](/helpers/genOtimesHouseholder.md): Generates two cells comprising of a set of complex Hermitian or Positive Semidefinite matrices of the respective dimensionalities.
- [genPos](/helpers/genPos.md): Generates a Positive Semidefinite matrix of the input dimensionality.