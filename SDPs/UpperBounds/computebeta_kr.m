function [rho_kr, beta_kr] = computebeta_kr(Pi, dim_B, k, PPT)
    cvx_begin sdp quiet
        variable rho_kr(size(Pi)) complex semidefinite;
        maximize HSIP(rho_kr, Pi)
        subject to
            trace(rho_kr) == 1; % density matrix constraint
            SymmetricExtension(rho_kr, k, dim_B, PPT, 0, eps^(1/4)); 
            % Symmetric extension constraint
            Tracenorm(Realignment(rho_kr)) <= 1; % realignment constraint
    cvx_end
    beta_kr = HSIP(rho_k, Pi);
end