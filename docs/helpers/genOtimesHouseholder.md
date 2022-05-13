# genOtimesHouseholder
Generates two cells comprising of a set of complex Hermitian or positive semidefinite matrices of the respective input dimensionalities

## Syntax
``genOtimesHouseholder(FILENAME, DIM_A, DIM_B, NUM_SUMMANDS, OPT)``

## Argument description
- ``FILENAME``: file location where the cells are stored.
- ``DIM_A``: integer that describes the dimension of the matrix for Alice.
- ``DIM_B``: integer that describes the dimension of the matrix for Bob.
- ``NUM_SUMMANDS``: number of matrices for Alice and Bob each
- ``OPT``: ``"Pos"`` if positive semidefinite matrix is required
       ``"Herm"`` if Hermitian matrix is required

## Example
    >> genOtimes("test.mat", 2, 2, 2, "Herm")

## Notes
The intention behind this function is to construct our desired matrix $\Pi$ as
{% raw %} 
$$
    \Pi = \sum\limits_{i=1}^{NUM_SUMMANDS} K_i \otimes L_i
$$
{% endraw %}
where
{% raw %} 
$
    \{K_i\}_{i=1}^{NUM\_SUMMANDS}
$
{% endraw %} 
is the set of
{% raw %} 
$
    DIM\_A \times DIM\_A
$
{% endraw %}
Hermitian or positive semidefinite matrices, and 
{% raw %}
$
    \{L_i\}_{i=1}^{NUM\_SUMMANDS}
$
{% endraw %} 
is the set of 
{% raw %}
$
    DIM\_B \times DIM\_B
$
{% endraw %}
Hermitian or positive semidefinite matrices.
- ``FILENAME`` also contains the time required to generate each of the ``NUM_SUMMANDS`` matrices.

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/helpers/genOtimesHouseholder.m)