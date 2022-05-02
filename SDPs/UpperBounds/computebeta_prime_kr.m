function [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr(Pi, dim_B, k, ...
                                                            PPT)
    cvx_begin sdp quiet
        variable rho_prime_kr(size(Pi)) complex semidefinite;
        maximize HSIP(rho_beta_prime_kr, Pi)
        subject to
            trace(rho_prime_kr) == 1; % density matrix
            SymmetricExtension(rho_prime_kr, k, dim_B, PPT, 1, eps^(1/4)); 
            % Bosonic symmetric extension constraint
            Tracenorm(Realignment(rho_prime_kr)) <= 1; 
            % realignment criterion
    cvx_end
    beta_prime_kr = HSIP(rho_beta_prime_kr, Pi);
end