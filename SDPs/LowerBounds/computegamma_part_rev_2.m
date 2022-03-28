function [sigma_A, sigma_B, gamma, i] = computegamma_part_rev_2(Pi, ...
                                            sigma, dim, threshold, ...
                                            N_seeSaw, opt)
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
    
    old_gamma = -Inf;
    for i = 1:N_seeSaw
        if strcmp(opt, "sigma_A")
            % Precompute sigma_A and find its P. e-vec
            C_B = PartialTrace(Tensor(sigma_A, eye(dim_B)) * Pi, 1);
            [V, ~] = eig(C_B);
            P_evec_C_B = V(:, end);
            sigma_B = P_evec_C_B * P_evec_C_B';
            sigma_B = sigma_B / trace(sigma_B); % normalize
            % Precompute sigma_B and find its P. e-vec
            C_A = PartialTrace(Tensor(eye(dim_A), sigma_B) * Pi, 2);
            [V, ~] = eig(C_A);
            P_evec_C_A = V(:, end);
            sigma_A = P_evec_C_A * P_evec_C_A';
            sigma_A = sigma_A / trace(sigma_A); % normalize
        else
            % Precompute sigma_B and find its P. e-vec
            C_A = PartialTrace(Tensor(eye(dim_A), sigma_B) * Pi, 2);
            [V, ~] = eig(C_A);
            P_evec_C_A = V(:, end);
            sigma_A = P_evec_C_A * P_evec_C_A';
            sigma_A = sigma_A / trace(sigma_A); % normalize
            % Precompute sigma_A and find its P. e-vec
            C_B = PartialTrace(Tensor(sigma_A, eye(dim_B)) * Pi, 1);
            [V, ~] = eig(C_B);
            P_evec_C_B = V(:, end);
            sigma_B = P_evec_C_B * P_evec_C_B';
            sigma_B = sigma_B / trace(sigma_B); % normalize
        end
            
        new_gamma = HSIP(Tensor(sigma_A, sigma_B), Pi);
        if abs(new_gamma - old_gamma) <= threshold
            break
        elseif new_gamma >= old_gamma
            old_gamma = new_gamma;
        end
    end
    gamma = new_gamma;
end