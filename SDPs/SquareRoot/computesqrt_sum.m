function [sigma_A, sigma_B, beta_sqrt, ...
          sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
          i] = computesqrt_sum(K_list, L_list, N_summands, dim_A, dim_B, ...
                                threshold, N_seeSaw, sigma_opt, ss_type)
                            
    cvx_begin sdp quiet
        variable W(dim_A * dim_B, dim_A * dim_B) complex semidefinite;
        variable sigma_A(dim_A, dim_A) complex semidefinite;
        variable sigma_B(dim_B, dim_B) complex semidefinite;
        % Objective function
        mu = 0;
        for j = 1:N_summands
            mu = mu + HSIP(W, Tensor(K_list{j}, L_list{j}));
        end
        maximize mu
        subject to
            [Tensor(sigma_A, eye(dim_A)) W; 
             W Tensor(eye(dim_B), sigma_B)] >= 0;
            trace(sigma_A) == 1; % density matrix constraint
            trace(sigma_B) == 1; % density matrix constraint
    cvx_end
    
    beta_sqrt = 0;
    for j = 1:N_summands
        prod = HSIP(sqrtm(sigma_A), K_list{j}) * HSIP(sqrt(sigma_B), L_list{j});
        beta_sqrt = beta_sqrt + prod;
    end
    
    % Polish our solution using see-saw
    if strcmp(sigma_opt, "sigma_A")
        if strcmp(ss_type, "simple")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_sum(K_list, L_list, N_summands, ...
                        sigma_A, dim_B, threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_1_sum(K_list, L_list, ...
                        N_summands, sigma_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_2_sum(K_list, L_list, ...
                        N_summands, sigma_A, dim_B, threshold, N_seeSaw, ...
                        "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", opt);
        end
    elseif strcmp(sigma_opt, "sigma_B")
        if strcmp(ss_type, "simple")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_sum(K_list, L_list, N_summands, ...
                        sigma_B, dim_A, threshold, N_seeSaw, "sigma_B");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_1_sum(K_list, L_list, ...
                        N_summands, sigma_B, dim_A, threshold, N_seeSaw, ...
                        "sigma_B");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_2_sum(K_list, L_list, ...
                        N_summands, sigma_B, dim_A, threshold, N_seeSaw, ...
                        "sigma_B");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", opt);
        end
    else
        error("Expected 'sigma_A' or 'sigma_B' for `sigma_opt`, got %s instead.", sigma_opt);
    end
end