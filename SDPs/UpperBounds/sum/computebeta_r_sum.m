%%  COMPUTEBETA_R_SUM   Computes the density matrix and the optimal value of the SDP involving the Realignment relaxation
%   The function has three required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices that have to be summed over
%
%   [RHO_R, BETA_R] = computebeta_r_sum(K_LIST, L_LIST, N_SUMMANDS)
%   computes the density matrix RHO_R and the optimal value BETA_R of the
%   SDP
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/sum/computebeta_r_sum.html
%
%   requires: cvx (http://cvxr.com/cvx/), Tensor
%   (http://qetlab.com/Tensor), TraceNorm (http://qetlab.com/TraceNorm),
%   Realignment (http://qetlab.com/Realignment), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_r, beta_r] = computebeta_r_sum(K_list, L_list, N_summands)
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_r(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_r, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            % density matrix constraint
            trace(rho_r) == 1;
            % realignment constraint
            TraceNorm(Realignment(rho_r)) <= 1;
    cvx_end
    beta_r = 0;
    for j = 1:N_summands
        beta_r = beta_r + HSIP(rho_r, Tensor(K_list{j}, L_list{j}));
    end
end