function [rho_PPT, beta_PPT] = computebeta_PPT_sum(K_list, L_list, ...
                                                    N_summands)
    dim = size(K_list{1}, 1) * size(L_list{1}, 1);
    cvx_begin sdp quiet
        variable rho_PPT(dim, dim) complex semidefinite;
        % Objective function
        obj = 0;
        for j = 1:N_summands
            obj = obj + HSIP(rho_PPT, Tensor(K_list{j}, L_list{j}));
        end
        maximize obj
        subject to
            trace(rho_PPT) == 1; % density matrix constraint
            PartialTranspose(rho_PPT, 2) >= 0; % PPT criterion 
    cvx_end
    beta_PPT = 0;
    for j = 1:N_summands
        beta_PPT = beta_PPT + HSIP(rho_PPT, Tensor(K_list{j}, L_list{j}));
    end
end