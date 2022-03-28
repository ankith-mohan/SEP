function [rho_k, beta_k] = computebeta_k_sum(K_list, L_list, N_summands, ...
                                                dim_B, k)
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_k(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_k, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            trace(rho_k) == 1; % density matrix constraint
            SymmetricExtension(rho_k, k, dim_B, 1, 0, eps^(1/4)); 
            % Symmetric extension constraint
    cvx_end
    beta_k = 0;
    for j = 1:N_summands
        beta_k = beta_k + HSIP(rho_k, Tensor(K_list{j}, L_list{j}));
    end
end