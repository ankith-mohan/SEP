function [rho_prime_k, beta_prime_k] = computebeta_prime_k(Pi, dim_B, k, ...
                                                            PPT)
    cvx_begin sdp quiet
        variable rho_prime_k(size(Pi)) complex semidefinite;
        maximize HSIP(rho_beta_prime_k, Pi)
        subject to
            trace(rho_prime_k) == 1; % density matrix
            SymmetricExtension(rho_prime_k, k, dim_B, PPT, 1, eps^(1/4)); 
            % Bosonic symmetric extension constraint
    cvx_end
    beta_prime_k = HSIP(rho_beta_prime_k, Pi);
end