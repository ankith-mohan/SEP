function [rho_r, beta_r] = computebeta_r_sum(K_list, L_list, N_summands)
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_r(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_r, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            trace(rho_r) == 1; % density matrix constraint
            TraceNorm(Realignment(rho_r)) <= 1; % realignment constraint 
    cvx_end
    beta_r = 0;
    for j = 1:N_summands
        beta_r = beta_r + HSIP(rho_r, Tensor(K_list{j}, L_list{j}));
    end
end