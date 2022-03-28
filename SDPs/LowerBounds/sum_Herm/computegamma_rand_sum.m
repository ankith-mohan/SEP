function [sigma_A_rand_best, sigma_B_rand_best, gamma_rand_best, ...
            i_best] = computegamma_rand_sum(K_list, L_list, N_summands, ...
                        dim_A, dim_B, threshold, N_seeSaw, N_rand, ss_type)
    gamma_rand_best = -Inf;
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
        
        if gamma_rand > gamma_rand_best
            sigma_A_rand_best = sigma_A_rand;
            sigma_B_rand_best = sigma_B_rand;
            gamma_rand_best = gamma_rand;
            i_best = i;
        end
    end
end