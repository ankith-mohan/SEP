%%  COMPUTEBETA_KR_SUM  Computes the density matrix and the optimal value for the SDP involving the Symmetric Extensions as well as the Realignment relaxations
%   This function has six required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_list: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices that have to be summed over
%       DIM_B: integer that describes the dimension of the matrix for Bob


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
            trace(rho_kr) == 1; % density matrix constraint
            SymmetricExtension(rho_kr, k, dim_B, PPT, 0, eps^(1/4)); 
            % Symmetric extension constraint
            TraceNorm(Realignment(rho_kr)) <= 1; % realignment constraint
    cvx_end
    beta_kr = 0;
    for j = 1:N_summands
        beta_kr = beta_kr + HSIP(rho_kr, Tensor(K_list{j}, L_list{j}));
    end
end