function [sigma_A, sigma_B, gamma, i] = computegamma_part(Pi, sigma, ...
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
            % Fix sigma_A and optimize sigma_B
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                maximize HSIP(Tensor(sigma_A, sigma_B), Pi)
                subject to
                    trace(sigma_B) == 1; % density matrix constraint
            cvx_end
            
            % Fix sigma_B and optimize sigma_A
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                maximize HSIP(Tensor(sigma_A, sigma_B), Pi)
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
        else
            %% Start from sigma_B
            % Fix sigma_B and optimize sigma_A
            cvx_begin sdp quiet
                variable sigma_A(dim_A, dim_A) complex semidefinite;
                maximize HSIP(Tensor(sigma_A, sigma_B), Pi)
                subject to
                    trace(sigma_A) == 1; % density matrix constraint
            cvx_end
            
            % Fix sigma_A and optimize sigma_B
            cvx_begin sdp quiet
                variable sigma_B(dim_B, dim_B) complex semidefinite;
                maximize HSIP(Tensor(sigma_A, sigma_B), Pi)
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