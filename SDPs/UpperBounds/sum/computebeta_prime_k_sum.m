%%  COMPUTEBETA_PRIME_K_SUM Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions relaxation
%   The function has six required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices to be summed over
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       K: integer that describes the number of symmetric extensions
%       PPT: 1 to impose PPT criterion
%            0 otherwise
%
%   [RHO_PRIME_K, BETA_PRIME_K] = computebeta_prime_k_sum(K_LIST, L_LIST, 
%   N_SUMMANDS, DIM_B, K, PPT) computes the density matrix RHO_PRIME_K and
%   the optimal value BETA_PRIME_K of the SDP
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/sum/computebeta_prime_k_sum.html
%
%   requires: cvx (http://cvxr.com/cvx/), Tensor
%   (http://qetlab.com/Tensor), SymmetricExtension
%   (http://qetlab.com/SymmetricExtension), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_prime_k, beta_prime_k] = computebeta_prime_k_sum(K_list, ...
                                            L_list, N_summands, dim_B, k, ...
                                            PPT)
    % Bosonic Extensions with or without PPT
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_prime_k(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_prime_k, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            % density matrix
            trace(rho_prime_k) == 1;
            % Bosonic symmetric extension constraint
            SymmetricExtension(rho_prime_k, k, dim_B, PPT, 1, eps^(1/4));
    cvx_end
    beta_prime_k = 0;
    for j = 1:N_summands
        prod = HSIP(rho_prime_k, Tensor(K_list{j}, L_list{j}));
        beta_prime_k = beta_prime_k + prod;
    end
end