function [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr_sum(K_list, ...
                                            L_list, N_summands, dim_B, k, ...
                                            PPT)
    % Bosonic Extensions with realignment and with or without PPT
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_prime_kr(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_prime_kr, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            trace(rho_prime_kr) == 1; % density matrix
            SymmetricExtension(rho_prime_kr, k, dim_B, PPT, 1, eps^(1/4)); 
            % Bosonic symmetric extension constraint
            TraceNorm(Realignment(rho_prime_kr)) <= 1; 
            % realignment constraint
    cvx_end
    beta_prime_kr = 0;
    for j = 1:N_summands
        prod = HSIP(rho_prime_kr, Tensor(K_list{j}, L_list{j}));
        beta_prime_kr = beta_prime_kr + prod;
    end
end