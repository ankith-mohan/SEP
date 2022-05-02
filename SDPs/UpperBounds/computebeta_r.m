function [rho_r, beta_r] = computebeta_r(Pi)
    cvx_begin sdp quiet
        variable rho_r(size(Pi)) complex semidefinite;
        maximize HSIP(rho_r, Pi)
        subject to
            trace(rho_r) == 1; % density matrix constraint
            TraceNorm(Realignment(rho_r)) <= 1; % realignment constraint 
    cvx_end
    beta_r = HSIP(rho_r, Pi);
end