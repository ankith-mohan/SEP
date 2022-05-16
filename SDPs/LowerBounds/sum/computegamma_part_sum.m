%%  COMPUTEGAMMA_PART_SUM   Computes the density matrix of the respective subsystems, the optimal value and the number of iterations for the see-saw algorithm
%   This function has eight required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices that have to be summed over
%       SIGMA: starting matrix for the see-saw algorithm
%       DIM: dimension of the starting matrix
%       THRESHOLD: threshold to determine convergence
%       N_SEESAW: maximum number of see-saw steps
%       OPT: "sigma_A" if starting matrix is for Alice's subsystem
%            "sigma_B" if starting matrix is for Bob's subsystem
%
%   [SIGMA_A, SIGMA_B, GAMMA, I] = computegamma_part_sum(K_LIST, L_LIST, 
%   N_SUMMANDS, SIGMA, DIM, THRESHOLD, N_SEESAW, OPT) computes the density
%   matrix of the respective subsystems SIGMA_A and SIGMA_B, the optimal
%   value GAMMA, and the number of see-saw steps before convergence I
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/LowerBounds/computegamma_part_sum.html
%
%   requires: cvx (http://cvxr.com/cvx), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [sigma_A, sigma_B, gamma, ...
            i] = computegamma_part_sum(K_list, L_list, N_summands, sigma, ...
                                        dim, threshold, N_seeSaw, opt)
    %% Initialization
    if strcmp(opt, "sigma_A")
        sigma_A = sigma;
        dim_A = size(sigma_A, 1);
        dim_B = dim;
    elseif strcmp(opt, "sigma_B")
        sigma_B = sigma;
        dim_A = dim;
        dim_B = size(sigma_B, 1);
    else
        error("Expected 'sigma_A' or 'sigma_B' for `opt`, got %s instead.", ...
              opt);
    end
    
    %% See-saw
    old_gamma = -Inf;
    for i = 1:N_seeSaw
        if strcmp(opt, "sigma_A")
            %% Start from sigma_A
            % Precompute <sigma_A, K_j> forall j and optimize sigma_B
            HSIP_A = zeros(1, N_summands);
            for j = 1:N_summands
                HSIP_A(j) = HSIP(sigma_A, K_list{j});
            end
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                % Objective function
                obj = 0;
                for j = 1:N_summands
                    obj = obj + HSIP_A(j) * HSIP(sigma_B, L_list{j});
                end
                maximize obj
                subject to
                    trace(sigma_B) == 1; % density matrix constraint
            cvx_end
            
            % Precompute <sigma_B, L_j> forall j and optimize sigma_A
            HSIP_B = zeros(1, N_summands);
            for j = 1:N_summands
                HSIP_B(j) = HSIP(sigma_B, L_list{j});
            end
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                % Objective function
                obj = 0;
                for j = 1:N_summands
                    obj = obj + HSIP(sigma_A, K_list{j}) * HSIP_B(j);
                end
                maximize obj
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
        else
            %% Start from sigma_B
            % Precompute <sigma_B, L_j> forall j and optimize sigma_A
            HSIP_B = zeros(1, N_summands);
            for j = 1:N_summands
                HSIP_B(j) = HSIP(sigma_B, L_list{j});
            end
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                % Objective function
                obj = 0;
                for j = 1:N_summands
                    obj = obj + HSIP(sigma_A, K_list{j}) * HSIP_B(j);
                end
                maximize obj
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
            
            % Precompute <sigma_A, K_j> forall j and optimize sigma_B
            HSIP_A = zeros(1, N_summands);
            for j = 1:N_summands
                HSIP_A(j) = HSIP(sigma_A, K_list{j});
            end
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                % Objective function
                obj = 0;
                for j = 1:N_summands
                    obj = obj + HSIP_A(j) * HSIP(sigma_B, L_list{j});
                end
                maximize obj
                subject to
                    trace(sigma_B) == 1; % density matrix constraint
            cvx_end
        end
        
        %% Update gamma
        new_gamma = 0;
        for j = 1:N_summands
            prod = HSIP(sigma_A, K_list{j}) * HSIP(sigma_B, L_list{j});
            new_gamma = new_gamma + prod;
        end
        
        %% Verify convergence
        if abs(new_gamma - old_gamma) <= threshold
            break
        elseif new_gamma > old_gamma
            old_gamma = new_gamma;
        end
    end
    gamma = new_gamma;
end