%% COMPUTEBETA_K    Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions relaxation
%   This function has four required input arguments:
%       PI: input matrix
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       K: integer that describes the number of symmetric extensions
%       PPT: 1 to impose PPT criterion
%            0 otherwise
%
%   [RHO_K, BETA_K] = computebeta_k(PI, DIM_B, K, PPT) computes the density
%   matrix RHO_K and the optimal value of the SDP BETA_K
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/computebeta_k.html
%
%   requires: cvx (http://cvxr.com/cvx/), SymmetricExtension (http://www.qetlab.com/SymmetricExtension), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_k, beta_k] = computebeta_k(Pi, dim_B, k, PPT)
    cvx_begin sdp quiet
        variable rho_k(size(Pi)) complex semidefinite;
        maximize HSIP(rho_k, Pi)
        subject to
            % density matrix constraint
            trace(rho_k) == 1;
            % Symmetric extension constraint
            SymmetricExtension(rho_k, k, dim_B, PPT, 0, eps^(1/4));
    cvx_end
    beta_k = HSIP(rho_k, Pi);
end