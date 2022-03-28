function [rho_k, beta_k] = computebeta_k(Pi, dim_B, k)
    cvx_begin sdp quiet
        variable rho_k(size(Pi)) complex semidefinite;
        maximize HSIP(rho_k, Pi)
        subject to
            trace(rho_k) == 1; % density matrix constraint
            SymmetricExtension(rho_k, k, dim_B, 1, 0, eps^(1/4)); 
            % Symmetric extension constraint
    cvx_end
    beta_k = HSIP(rho_k, Pi);
end