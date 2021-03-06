%% GENPOS   Generates a complex positive semidefinite matrix of the input dimensionality
%   This function has one required input argument:
%       DIM: integer that describes the local dimension of the matrix
%
%   P = genPos(DIM) is a positive semidefinite matrix of size DIM x DIM
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/genPos.html
%
%   requires: IsPSD (http://www.qetlab.com/IsPSD)
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function P = genPos(dim)
    % Create a Hermitian matrix
    H = genHerm(dim);
    % Make it PSD
    P = H' * H;
    % Sure its PSD?
    assert(IsPSD(P), "P is not PSD");
end