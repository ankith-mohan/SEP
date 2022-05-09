%% HSIP Computes the Hilbert-Schmidt inner product between the input matrices
%   This function requires two input arguments:
%       A: matrix of size m x n
%       B: matrix of size m x p
%
%   HSIP(A, B) is the Hilbert-Schmidt inner product of A and B
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/HSIP.html
%
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function dot = HSIP(A, B)
    dot = real(trace(A' * B));
end