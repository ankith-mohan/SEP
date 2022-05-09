%%  GENHERM  Generates a complex Hermitian matrix of the input dimensionality
%   This function has one required input argument:
%       DIM: integer that describes the local dimension of the matrix
%
%   H = genHerm(DIM) is a Hermitian matrix of size dim x dim
%
%   These matrices can be used in CVX or QETLAB.
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/genHerm.html
%
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function H = genHerm(dim)
    % Create a complex matrix
    M = complex(rand(dim), rand(dim));
    % Make it Hermitian
    H = 0.5 * (M + M');
    % Are you really sure its Hermitian?
    assert(ishermitian(H), "H is not Hermitian");
end