%%  COMPUTEBETA_PPT_SUM  Computes the density matrix and the optimal value for the SDP involving the PPT relaxation
%   The function has three required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices that have to be summed over
%
%   [RHO_PPT, BETA_PPT] = computebeta_PPT_sum(K_LIST, L_LIST, N_SUMMANDS)
%   computes the density matrix RHO_PPT and the optimal value BETA_PPT of
%   the SDP
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/sum/computebeta_PPT_sum.html
%
%   requires: cvx (http://cvxr.com/cvx/), Tensor
%   (http://qetlab.com/Tensor), PartialTranspose
%   (http://qetab.com/PartialTranspose), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_PPT, beta_PPT] = computebeta_PPT_sum(K_list, L_list, ...
                                                    N_summands)
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_PPT(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_PPT, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            % density matrix constraint
            trace(rho_PPT) == 1;
            % PPT criterion 
            PartialTranspose(rho_PPT, 2) >= 0;
    cvx_end
    beta_PPT = 0;
    for j = 1:N_summands
        beta_PPT = beta_PPT + HSIP(rho_PPT, Tensor(K_list{j}, L_list{j}));
    end
end