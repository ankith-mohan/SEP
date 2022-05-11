%%  COMPUTEBETA_PRIME_K Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions 
%   This function has four required input arguments:
%       PI: input matrix
%       DIM_B: input that describes the dimension of the matrix for Bob
%       K: integer that describes the number of bosonic extensions
%       PPT: 1 to impose PPT criterion
%            0 otherwise
%
%   [RHO_PRIME_K, BETA_PRIME_K] = computebeta_prime_k(PI, DIM_B, K, PPT)
%   computes the density matrix RHO_PRIME_K and the optimal value of the
%   SDP BETA_PRIME_K
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/computebeta_prime_k.html
%
%   requires: cvx (http://cvxr.com/cvx), SymmetricExtension
%   (http://www.qetlab.com/SymmetricExtension), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_prime_k, beta_prime_k] = computebeta_prime_k(Pi, dim_B, k, ...
                                                            PPT)
    cvx_begin sdp quiet
        variable rho_prime_k(size(Pi)) complex semidefinite;
        maximize HSIP(rho_prime_k, Pi)
        subject to
            % density matrix
            trace(rho_prime_k) == 1;
            % Bosonic symmetric extension constraint
            SymmetricExtension(rho_prime_k, k, dim_B, PPT, 1, eps^(1/4));
    cvx_end
    beta_prime_k = HSIP(rho_prime_k, Pi);
end