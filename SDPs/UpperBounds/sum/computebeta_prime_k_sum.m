function [rho_prime_k, beta_prime_k] = computebeta_prime_k_sum(K_list, ...
                                            L_list, N_summands, dim_B, k, ...
                                            PPT)
    % Bosonic Extensions with or without PPT
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_prime_k(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_prime_k, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            trace(rho_prime_k) == 1; % density matrix
            SymmetricExtension(rho_prime_k, k, dim_B, PPT, 1, eps^(1/4)); 
            % Bosonic symmetric extension constraint
    cvx_end
    beta_prime_k = 0;
    for j = 1:N_summands
        prod = HSIP(rho_prime_k, Tensor(K_list{j}, L_list{j}));
        beta_prime_k = beta_prime_k + prod;
    end
end