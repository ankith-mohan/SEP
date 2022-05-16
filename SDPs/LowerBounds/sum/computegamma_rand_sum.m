%%  COMPUTEGAMMA_RAND_SUM   Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points
%   This function has ten required input arguments:
%       K_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_A x DIM_A
%       L_LIST: cell of length N_SUMMANDS where each entry is a matrix of
%       size DIM_B x DIM_B
%       N_SUMMANDS: number of matrices that have to be summed over
%       DIM_A: dimension of Alice's subsystem
%       DIM_B: dimension of Bob's subsystem
%       THRESHOLD: threshold to determine convergence
%       N_SEESAW: maximum number of see-saw steps
%       N_RAND: number of random starting points
%       SS_TYPE: "simple" if COMPUTEGAMMA_PART
%                "rev_1" if COMPUTEGAMMA_PART_REV_1
%                "rev_2" if COMPUTEGAMMA_PART_REV_2
%       PLUS_TWO: if both the maximally mixed state and the state of
%       uniform superposition must be added to the list of starting points
%
%   [SIGMA_A_RAND_BEST, SIGMA_B_RAND_BEST, GAMMA_RAND_BEST, I_BEST, 
%   GAMMA_RAND_LIST] = computegamma_rand(K_LIST, L_LIST, N_SUMMANDS, DIM_A, 
%   DIM_B, THRESHOLD, N_SEESAW, N_RAND, SS_TYPE, PLUS_TWO) returns the 
%   density matrix of the respective subsystems SIGMA_A_RAND_BEST and 
%   SIGMA_B_RAND_BEST, corresponding to the best optimal value 
%   GAMMA_RAND_BEST, the number of see-saw steps before convergence I_BEST, 
%   and the optimal values for all the starting points GAMMA_RAND_LIST
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/LowerBounds/computegamma_rand_sum.html
%
%   requires: computegamma_part_sum.m, computegamma_part_rev_1_sum.m, 
%   computegamma_part_rev_2_sum.m, genHouseholder.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, ...
            i_best, gamma_rand_list] = computegamma_rand_sum(K_list, L_list, ...
                                        N_summands, dim_A, dim_B, threshold, ...
                                        N_seeSaw, N_rand, ss_type, plus_two)
    gamma_rand_best = -Inf;
    gamma_rand_list = zeros(1, N_rand+plus_two);
    
    if plus_two == 2
        % Maximally-mixed state
        if strcmp(ss_type, "simple")
            [sigma_A_MM, sigma_B_MM, gamma_MM, ...
                i] = computegamma_part_sum(K_list, L_list, N_summands, ...
                        eye(dim_A)/dim_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_MM, sigma_B_MM, gamma_MM, ...
                i] = computegamma_part_rev_1_sum(K_list, L_list, N_summands, ...
                        eye(dim_A)/dim_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_MM, sigma_B_MM, gamma_MM, ...
                i] = computegamma_part_rev_2_sum(K_list, L_list, N_summands, ...
                        eye(dim_A)/dim_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", ss_type);
        end
        gamma_rand_list(1) = gamma_MM;
        if gamma_MM > gamma_rand_best
            sigma_A_rand_best = sigma_A_MM;
            sigma_B_rand_best = sigma_B_MM;
            gamma_rand_best = gamma_MM;
            i_best = i;
        end

        % Uniform superposition state
        if strcmp(ss_type, "simple")
            [sigma_A_US, sigma_B_US, gamma_US, ...
                i] = computegamma_part_sum(K_list, L_list, N_summands, ...
                        ones(dim_A)/dim_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_US, sigma_B_US, gamma_US, ...
                i] = computegamma_part_rev_1_sum(K_list, L_list, N_summands, ...
                        ones(dim_A)/dim_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_US, sigma_B_US, gamma_US, ...
                i] = computegamma_part_rev_2_sum(K_list, L_list, N_summands, ...
                        ones(dim_A)/dim_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", ss_type);
        end
        gamma_rand_list(2) = gamma_US;
        if gamma_US > gamma_rand_best
            sigma_A_rand_best = sigma_A_US;
            sigma_B_rand_best = sigma_B_US;
            gamma_rand_best = gamma_US;
            i_best = i;
        end
    end

    for r = 1:N_rand
        % sigma_A_rand is random PSD normalized by its trace norm
        sigma_A_rand = genHouseholder(dim_A, "Pos"); % PSD
        sigma_A_rand = sigma_A_rand / trace(sigma_A_rand); % normalized PSD
        
        if strcmp(ss_type, "simple")
            [sigma_A_rand, sigma_B_rand, gamma_rand, ...
                i] = computegamma_part_sum(K_list, L_list, N_summands, ...
                        sigma_A_rand, dim_B, threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_rand, sigma_B_rand, gamma_rand, ...
                i] = computegamma_part_rev_1_sum(K_list, L_list, ...
                        N_summands, sigma_A_rand, dim_B, threshold, ...
                        N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_rand, sigma_B_rand, gamma_rand, ...
                i] = computegamma_part_rev_2_sum(K_list, L_list, ...
                        N_summands, sigma_A_rand, dim_B, threshold, ...
                        N_seeSaw, "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", ss_type);
        end
        gamma_rand_list(r+plus_two) = gamma_US;
        if gamma_rand > gamma_rand_best
            sigma_A_rand_best = sigma_A_rand;
            sigma_B_rand_best = sigma_B_rand;
            gamma_rand_best = gamma_rand;
            i_best = i;
        end
    end
end