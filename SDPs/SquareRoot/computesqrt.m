function [sigma_A, sigma_B, beta_sqrt, ...
          sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
          i] = computesqrt(Pi, dim_A, dim_B, threshold, N_seeSaw, ...
                           sigma_opt, ss_type)
    cvx_begin sdp quiet
        variable W(size(Pi)) complex semidefinite;
        variable sigma_A(dim_A, dim_A) complex semidefinite;
        variable sigma_B(dim_B, dim_B) complex semidefinite;
        maximize HSIP(W, Pi) % mu
        subject to
            [Tensor(sigma_A, eye(dim_A)) W; 
             W Tensor(eye(dim_B), sigma_B)] >= 0;
            trace(sigma_A) == 1; % density matrix constraint
            trace(sigma_B) == 1; % density matrix constraint
    cvx_end
    
    beta_sqrt = HSIP(Tensor(sqrtm(sigma_A), sqrtm(sigma_B)), Pi);
    
    % Polish our solution using see-saw
    if strcmp(sigma_opt, "sigma_A")
        if strcmp(ss_type, "simple")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part(Pi, sigma_A, dim_B, threshold, ...
                                        N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_1(Pi, sigma_A, dim_B, ...
                                        threshold, N_seeSaw, "sigma_A");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_2(Pi, sigma_A, dim_B, ...
                                        threshold, N_seeSaw, "sigma_A");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", opt);
        end
    elseif strcmp(sigma_opt, "sigma_B")
        if strcmp(ss_type, "simple")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part(Pi, sigma_B, dim_A, threshold, ...
                                        N_seeSaw, "sigma_B");
        elseif strcmp(ss_type, "rev_1")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_1(Pi, sigma_B, dim_A, ...
                                        threshold, N_seeSaw, "sigma_B");
        elseif strcmp(ss_type, "rev_2")
            [sigma_A_sqrt, sigma_B_sqrt, gamma_sqrt, ...
                i] = computegamma_part_rev_2(Pi, sigma_B, dim_A, ...
                                        threshold, N_seeSaw, "sigma_B");
        else
            error("Expected 'simple', 'rev_1' or 'rev_2' for `ss_type`, got %s instead.", opt);
        end
    else
        error("Expected 'sigma_A' or 'sigma_B' for `sigma_opt`, got %s instead.", sigma_opt);
    end
end