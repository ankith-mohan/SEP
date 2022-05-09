# genHouseholder
Generates a complex Hermitian or Positive Semidefinite matrix of the input dimensionality

## Syntax
``H = genHouseholder(DIM, OPT)``

## Argument description
- ``DIM``: integer that describes the local dimension of the matrix.
- ``OPT``: ``"Pos"`` if positive semidefinite matrix is required
       ``"Herm"`` if Hermitian matrix is required

## Example
    >> genHouseholder(2, "Pos")

    ans = 

        0.5621 + 0.0000i    -0.0205 + 0.0047i
        -0.0205 - 0.0047i   0.5107 + 0.0000i

    >> genHouseholder(2, "Herm")

    ans = 

        0.7157 + 0.0000i    -0.1227 + 0.1079i
        -0.1227 - 0.1079i   0.7877 + 0.0000i

## Notes
This function does the following:
1. The spectrum is ``DIM`` numbers sampled uniformly in the range: [0, 1] if ``OPT`` is ``"Pos"``, [-1, 1] if ``OPT`` is ``"Herm"``
2. This spectrum forms the diagonal entries of the matrix $$\Lambda$$
3. Creates a random complex matrix ``M``
4. Orthogonal ``M`` to get the matrix ``Q``
5. ``M = Q *`` $$\Lambda$$ ``* Q``
6. ``M = 0.5 * (M + M')`` to make ``M`` Hermitian
7. ``M = M' * M`` if ``OPT`` is ``"Pos"``

## Source code
[Source](https://github.com/ankith-mohan/SEP/blob/main/helpers/genHouseholder.m)