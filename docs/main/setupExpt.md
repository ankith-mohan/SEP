# setupExpt
Computes the upper and lower bounds for ``N_OPS`` instances

## Syntax
``setupExpt(DIM_A, DIM_B, K_SE, K_BSE, N_SEESAW, N_RAND, SUM, NUM_SUMMANDS, N_OPS, THRESHOLD, OPT, FILENAME, PROJ, DONT_DO, SHOULD_SAVE, PLUS_TWO_RAND)``

## Argument description
#### Input arguments
- ``DIM_A``: integer that describes the dimension of the matrix for Alice
- ``DIM_B``: integer that describes the dimension of the matrix for Bob
- ``K_SE``: integer that describes the maximum number of symmetric extensions
- ``K_BSE``: integer that describes the maximum number of bosonic extensions
- ``N_SEESAW``: maximum number of see-saw steps
- ``N_RAND``: number of random starting points
- ``SUM``:
    - 1 if ``GENOTIMESHOUSEHOLDER`` is used to generate the cells of matrices
    - 0 otherwise
- ``NUM_SUMMANDS``: number of matrices that have to be summed over
- ``N_OPS``: number of instances of the input matrix to compute the bounds for
- ``THRESHOLD``: threshold to determine convergence
- ``OPT``:
    - ``"sigma_A"`` if starting matrix is for Alice's subsystem
    - ``"sigma_B"`` if starting matrix is for Bob's subsystem
- ``FILENAME``: ``".mat"`` file where all the results have to be written
- ``PROJ``:
    - 0 to use ``GENHOUSEHOLDER``
    - 1 to use ``GENHERM`` if ``OPT == "Herm"``
               ``GENPOS`` if ``OPT == "Pos"``
- ``DONT_DO``: list that contains what functions must NOT be called
    - ``"PPT"``: ``COMPUTEBETA_PPT`` if ``SUM == 0`` else ``COMPUTEBETA_PPT_SUM``
    - ``"realignment"``: ``COMPUTEBETA_R`` if ``SUM == 0`` else ``COMPUTEBETA_R_SUM``
    - ``"SymExtOnly"``: ``COMPUTEBETA_K(PPT = 0)`` if ``SUM == 0`` else ``COMPUTEBETA_K_SUM(PPT = 0)``
    - ``"SymExtOnlyrlg"``: ``COMPUTEBETA_KR(PPT = 0)`` if ``SUM == 0`` else ``COMPUTEBETA_KR_SUM(PPT = 0)``
    - ``"SymExt"``: ``COMPUTEBETA_K(PPT = 1)`` if ``SUM == 0`` else ``COMPUTEBETA_K_SUM(PPT = 1)``
    - ``"SymExtrlg"``: ``COMPUTEBETA_KR(PPT = 1)`` if ``SUM == 0`` else ``COMPUTEBETA_KR_SUM(PPT = 1)``
    - ``"BosExtOnly"``: ``COMPUTEBETA_PRIME_K(PPT = 0)`` if ``SUM == 0`` else ``COMPUTEBETA_PRIME_K_SUM(PPT = 0)``
    - ``"BosExtOnlyrlg"``: ``COMPUTEBETA_PRIME_KR(PPT = 0)`` if ``SUM == 0`` else ``COMPUTEBETA_PRIME_KR_SUM(PPT = 0)``
    - ``"BosExt"``: ``COMPUTEBETA_PRIME_K(PPT = 1)`` if ``SUM == 0`` else ``COMPUTEBETA_PRIME_K_SUM(PPT = 1)``
    - ``"BosExtrlg"``: ``COMPUTEBETA_PRIME_KR(PPT = 1)`` if ``SUM == 0`` else ``COMPUTEBETA_PRIME_KR_SUM(PPT = 1)``
    - ``"MM_simple"``: ``COMPUTEGAMMA_PART(SIGMA = eye(DIM_A)/DIM_A)`` if ``SUM == 0`` else ``COMPUTEGAMMA_PART_SUM(SIGMA = eye(DIM_A)/DIM_A)``
    - ``"MM_rev_1"``: ``COMPUTEGAMMA_PART_REV_1(SIGMA = eye(DIM_A)/DIM_A)`` if ``SUM == 0`` else ``COMPUTEGAMMA_PART_REV_1_SUM(SIGMA = eye(DIM_A)/DIM_A)``
    - ``"MM_rev_2"``: ``COMPUTEGAMMA_PART_REV_2(SIGMA = eye(DIM_A)/DIM_A)`` if ``SUM == 0`` else ``COMPUTEGAMMA_PART_REV_2_SUM(SIGMA = eye(DIM_A)/DIM_A)``
    - ``"US_simple"``: ``COMPUTEGAMMA_PART(SIGMA = ones(DIM_A)/DIM_A)`` if ``SUM == 0`` else ``COMPUTEGAMMA_PART_SUM(SIGMA = ones(DIM_A)/DIM_A)``
    - ``"US_rev_1"``: ``COMPUTEGAMMA_PART_REV_1(SIGMA = ones(DIM_A)/DIM_A)`` if ``SUM == 0`` else ``COMPUTEGAMMA_PART_REV_1_SUM(SIGMA = ones(DIM_A)/DIM_A)``
    - ``"US_rev_2"``: ``COMPUTEGAMMA_PART_REV_2(SIGMA = ones(DIM_A)/DIM_A)`` if ``SUM == 0`` else ``COMPUTEGAMMA_PART_REV_2_SUM(SIGMA = ones(DIM_A)/DIM_A)``
    - ``"rand_simple"``: ``COMPUTEGAMMA_RAND(SS_TYPE = "simple")`` if ``SUM == 0`` else ``COMPUTEGAMMA_RAND_SUM(SS_TYPE = "simple")``
    - ``"rand_rev_1"``: ``COMPUTEGAMMA_RAND(SS_TYPE = "rev_1")`` if ``SUM == 0`` else ``COMPUTEGAMMA_RAND_SUM(SS_TYPE = "rev_1")``
    - ``"MM_rev_2"``: ``COMPUTEGAMMA_RAND(SS_TYPE = "rev_2")`` if ``SUM == 0`` else ``COMPUTEGAMMA_RAND_SUM(SS_TYPE = "rev_2")``
- ``SHOULD_SAVE``:
    - 1 if the density matrices returned by each of the functions need to be stored in ``FILENAME``
    - 0 otherwise
- ``PLUS_TWO_RAND``:
    - 2 if both the maximally mixed state and the state of uniform superposition must be added to the list of starting points
    - 0 otherwise

## Example
    >> setupExpt(2, 2, 2, 2, 100, 100, 1, 10, 100, 1e-13, "Herm", "Alice=2,Bob=2,Herm", 0, [], 0, 2)

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/main/setupExpt.m)