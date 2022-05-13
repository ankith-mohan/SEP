%%  COMPUTEBETA_KR_SUM  Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations
%   This function has six required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_list: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices that have to be summed over
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       K: integer that describes the number of symmetric extensions
%       PPT: 1 to impose PPT criterion
%            0 otherwise
%
%   [RHO_KR, BETA_KR] = computebeta_kr_sum(K_LIST, L_LIST, N_SUMMANDS,
%   DIM_B, K, PPT) computes the density matrix RHO_KR and the optimal value of the SDP BETA_KR
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/sum/computebeta_kr_sum.html
%
%   requires: cvx (http://cvxr.com/cvx/), Tensor
%   (http://qetlab.com/Tensor), SymmetricExtension
%   (http://qetlab.com/SymmetricExtension), TraceNorm
%   (http://qetlab.com/TraceNorm), Realignment
%   (http://qetlab.com/Realignment), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_kr, beta_kr] = computebeta_kr_sum(K_list, L_list, N_summands, ...
                                                dim_B, k, PPT)
    % Symmetric Extensions with realignment and with or without PPT
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_kr(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_kr, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            % density matrix constraint
            trace(rho_kr) == 1;
            % Symmetric extension constraint
            SymmetricExtension(rho_kr, k, dim_B, PPT, 0, eps^(1/4));
            % realignment constraint
            TraceNorm(Realignment(rho_kr)) <= 1;
    cvx_end
    beta_kr = 0;
    for j = 1:N_summands
        beta_kr = beta_kr + HSIP(rho_kr, Tensor(K_list{j}, L_list{j}));
    end
end