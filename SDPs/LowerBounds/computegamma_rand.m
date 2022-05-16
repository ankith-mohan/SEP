%%  COMPUTEGAMMA_RAND   Computes the density matrix of the respective subsystems, largest of the optimal values, the number of iterations, and the list of optimal values for the see-saw algorithm with random starting points
%   This function has eight required input arguments:
%       PI: input matrix
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
%   GAMMA_RAND_LIST] = computegamma_rand(PI, DIM_A, DIM_B, THRESHOLD,
%   N_SEESAW, N_RAND, SS_TYPE, PLUS_TWO) returns the density matrix of the
%   respective subsystems SIGMA_A_RAND_BEST and SIGMA_B_RAND_BEST, 
%   corresponding to the best optimal value GAMMA_RAND_BEST, the number of 
%   see-saw steps before convergence I_BEST, and the optimal values for all
%   the starting points GAMMA_RAND_LIST
%
%   URL:
%   https://ankith-mohan.github.io/SEP/SDPs/LowerBounds/computegamma_rand.html
%
%   requires: computegamma_part.m, computegamma_part_rev_1.m, 
%   computegamma_part_rev_2.m, genHouseholder.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, ...
            i_best, gamma_rand_list] = computegamma_rand(Pi, dim_A, dim_B, ...
                                        threshold, N_seeSaw, N_rand, ...
                                        ss_type, plus_two)
    gamma_rand_best = -Inf;
    gamma_rand_list = [];

    if plus_two == 1
        % Maximally-mixed state
        if strcmp(ss_type, "simple")
            [sigma_A_MM, sigma_B_MM, gamma_MM, ...
                i] = computegamma_part(Pi, eye(dim_A)/dim_A, dim_B, ...
                                        threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_MM, sigma_B_MM, gamma_MM, ...
                i] = computegamma_part_rev_1(Pi, eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_MM, sigma_B_MM, gamma_MM, ...
                i] = computegamma_part_rev_2(Pi, eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", ss_type);
        end
        gamma_rand_list = [gamma_rand_list gamma_MM];
        if gamma_MM > gamma_rand_best
            sigma_A_rand_best = sigma_A_MM;
            sigma_B_rand_best = sigma_B_MM;
            gamma_rand_best = gamma_MM;
            i_best = i;
        end

        % Uniform superposition state
        if strcmp(ss_type, "simple")
            [sigma_A_US, sigma_B_US, gamma_US, ...
                i] = computegamma_part(Pi, ones(dim_A)/dim_A, dim_B, ...
                                        threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_US, sigma_B_US, gamma_US, ...
                i] = computegamma_part_rev_1(Pi, ones(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_US, sigma_B_US, gamma_US, ...
                i] = computegamma_part_rev_2(Pi, ones(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", ss_type);
        end
        gamma_rand_list = [gamma_rand_list gamma_US];
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
                i] = computegamma_part(Pi, sigma_A_rand, dim_B, threshold, ...
                                        N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_rand, sigma_B_rand, gamma_rand, ...
                i] = computegamma_part_rev_1(Pi, sigma_A_rand, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_rand, sigma_B_rand, gamma_rand, ...
                i] = computegamma_part_rev_2(Pi, sigma_A_rand, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", ss_type);
        end
        gamma_rand_list = [gamma_rand_list gamma_rand];
        if gamma_rand > gamma_rand_best
            sigma_A_rand_best = sigma_A_rand;
            sigma_B_rand_best = sigma_B_rand;
            gamma_rand_best = gamma_rand;
            i_best = i;
        end
    end
end