# genMat
Generates a complex Hermitian or positive semidefinite matrix of the input dimensionality normalized by its trace norm

## Syntax
``M = genMat(DIM, OPT, PROJ)``

## Argument description
- ``DIM``: integer that describes the local dimension of the matrix.
- ``OPT``: ``"Pos"`` if positive semidefinite matrix is required
           ``"Herm"`` if Hermitian matrix is required
- ``PROJ``: 0 to use ``genHouseholder``
            1 to use ``genHerm`` if ``OPT == "Herm"`` or ``genPos`` if ``OPT == "Pos"``

## Example
    >> genMat(2, "Pos", 0)

    ans =

        0.1487 + 0.0000i  -0.0248 - 0.0080i
        -0.0248 + 0.0080i   0.3185 + 0.0000i

    >> genMat(2, "Herm", 0)

    ans =

        -0.3591 + 0.0000i  -0.1420 - 0.0211i
        -0.1420 + 0.0211i  -0.9428 + 0.0000i

    >> genMat(2, "Pos", 1)

    ans = 

        0.7297 + 0.0000i   0.3453 - 0.2192i
        0.3453 + 0.2192i   0.3809 + 0.0000i
        
    >> genMat(2, "Herm", 1)

    ans =

        0.4443 + 0.0000i   0.4434 - 0.3951i
        0.4434 + 0.3951i   0.3654 + 0.0000i

## Source code
[Source](https://github.com/ankith-mohan/SEP/blob/main/helpers/genMat.m)