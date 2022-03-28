function [sigma_A, sigma_B, gamma, ...
            i] = computegamma_part_rev_1_sum(K_list, L_list, N_summands, ...
                                    sigma, dim, threshold, N_seeSaw, opt)
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
        error("Expected 'sigma_A' or 'sigma_B' for `opt`, got %s instead.", opt);
    end
    
    %% See-saw
    old_gamma = -Inf;
    for i = 1:N_seeSaw
        if strcmp(opt, "sigma_A")
            %% Start from sigma_A
            % Precompute sigma_A and optimize sigma_B
            C_B = zeros('like', L_list{1});
            for j = 1:N_summands
                C_B = C_B + HSIP(sigma_A, K_list{j}) * L_list{j};
            end
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                maximize HSIP(C_B, sigma_B)
                subject to
                    trace(sigma_B) == 1; % density matrix constraint
            cvx_end
            
            % Precompute sigma_B and optimize sigma_A
            C_A = zeros('like', K_list{1});
            for j = 1:N_summands
                C_A = C_A + HSIP(sigma_B, L_list{j}) * K_list{j};
            end
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                maximize HSIP(C_A, sigma_A)
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
        else
            %% Start from sigma_B
            % Precompute sigma_B and optimize sigma_A
            C_A = zeros('like', K_list{1});
            for j = 1:N_summands
                C_A = C_A + HSIP(sigma_B, L_list{j}) * K_list{j};
            end
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                maximize HSIP(C_A, sigma_A)
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
            
            % Precompute sigma_A and optimize sigma_B
            C_B = zeros('like', L_list{1});
            for j = 1:N_summands
                C_B = C_B + HSIP(sigma_A, K_list{j}) * L_list{j};
            end
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                maximize HSIP(C_B, sigma_B)
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
        elseif new_gamma >= old_gamma
            old_gamma = new_gamma;
        end
    end
    gamma = new_gamma;
end