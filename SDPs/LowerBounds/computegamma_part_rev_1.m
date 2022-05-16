%%  COMPUTEGAMMA_PART_REV_1 Computes the density matrix of the respective subsystems, optimal value and the number of iterations for revision 1 of the see-saw algorithm
%   This function has six required input arguments:
%       PI: input matrix
%       SIGMA: Starting matrix for the see-saw algorithm
%       DIM: dimension of the starting matrix
%       THRESHOLD: threshold to determine convergence
%       N_SEESAW: maximum number of see-saw steps
%       OPT: "sigma_A" if starting matrix is for Alice's subsystem
%            "sigma_B" if starting matrix is for Bob's subsystem
%
%   [SIGMA_A, SIGMA_B, GAMMA, I] = computegamma_part_rev_1(PI, SIGMA, DIM,
%   THRESHOLD, N_SEESAW, OPT) computes the density matrix of the
%   respective subsystems SIGMA_A and SIGMA_B, the optimal value GAMMA, and
%   the number of see-saw steps before convergence I
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/LowerBounds/computegamma_part_rev_1.html
%
%   requires: cvx (http://cvxr.com/cvx), PartialTrace (http://qetlab.com/PartialTrace), 
%   Tensor (http://qetlab.com/Tensor), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [sigma_A, sigma_B, gamma, i] = computegamma_part_rev_1(Pi, ...
                                            sigma, dim, threshold, ...
                                            N_seeSaw, opt)
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
            % Precompute sigma_A and optimize sigma_B
            C_B = PartialTrace(Tensor(sigma_A, eye(dim_B)) * Pi, 1);
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                maximize HSIP(C_B, sigma_B)
                subject to
                    trace(sigma_B) == 1; % density matrix constraint
            cvx_end
            
            % Precompute sigma_B and optimize sigma_A
            C_A = PartialTrace(Tensor(eye(dim_A), sigma_B) * Pi, 2);
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                maximize HSIP(sigma_A, C_A)
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
        else
            %% Start from sigma_B
            % Precompute sigma_B and optimize sigma_A
            C_A = PartialTrace(Tensor(eye(dim_A), sigma_B) * Pi, 2);
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                maximize HSIP(sigma_A, C_A)
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
            
            % Precompute sigma_A and optimize sigma_B
            C_B = PartialTrace(Tensor(sigma_A, eye(dim_B)) * Pi, 1);
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                maximize HSIP(C_B, sigma_B)
                subject to
                    trace(sigma_B) == 1; % density matrix constraint
            cvx_end
        end
        
        %% Update gamma
        new_gamma = HSIP(Tensor(sigma_A, sigma_B), Pi);
        
        %% Verify convergence
        if abs(new_gamma - old_gamma) <= threshold
            break
        elseif new_gamma >= old_gamma
            old_gamma = new_gamma;
        end
    end
    gamma = new_gamma;
end