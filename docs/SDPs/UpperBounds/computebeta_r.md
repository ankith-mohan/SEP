# computebeta_r
Computes the density matrix and the optimal value for the SDP involving the Realignment relaxation

## Syntax
``[RHO_R, BETA_R] = computebeta_r(PI)``

## Argument description
#### Input arguments
- ``PI``: input matrix

#### Output arguments
- ``RHO_R``: density matrix of the Realignment SDP
- ``BETA_R``: optimal value of the Realignment SDP

## Example
    >> [rho_r, beta_r] = computebeta_r(genMat(4, "Herm", 0))

    rho_r =

        0.2885 + 0.0000i   0.1517 + 0.0037i   0.3761 + 0.0363i   0.1973 + 0.0239i
        0.1517 - 0.0037i   0.0798 + 0.0000i   0.1982 + 0.0142i   0.1041 + 0.0100i
        0.3761 - 0.0363i   0.1982 - 0.0142i   0.4948 + 0.0000i   0.2602 + 0.0064i
        0.1973 - 0.0239i   0.1041 - 0.0100i   0.2602 - 0.0064i   0.1369 + 0.0000i

    beta_r =

        0.9635

    >> [rho_r, beta_r] = computebeta_r(genMat(4, "Pos", 0))

    rho_r =

        0.2300 + 0.0000i   0.1441 + 0.0160i   0.1524 - 0.2974i   0.1162 - 0.1757i
        0.1441 - 0.0160i   0.0914 + 0.0000i   0.0747 - 0.1970i   0.0606 - 0.1182i
        0.1524 + 0.2974i   0.0747 + 0.1970i   0.4856 + 0.0000i   0.3043 + 0.0339i
        0.1162 + 0.1757i   0.0606 + 0.1182i   0.3043 - 0.0339i   0.1930 + 0.0000i

    beta_r =

        0.7304

## Source Code
[Source](https://github.com/ankith-mohan/SEP/blob/main/helpers/computebeta_r.m)

## References
1. Chen, K., & Wu, L. A. (2002). A matrix realignment method for recognizing entanglement. arXiv preprint quant-ph/0205017.