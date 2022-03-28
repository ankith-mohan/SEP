function [rho_PPT, beta_PPT] = computebeta_PPT(Pi)
    cvx_begin sdp quiet
        variable rho_PPT(size(Pi)) complex semidefinite;
        maximize HSIP(rho_PPT, Pi)
        subject to
            trace(rho_PPT) == 1; % density matrix constraint
            PartialTranspose(rho_PPT, 2) >= 0; % PPT criterion 
    cvx_end
    beta_PPT = HSIP(rho_PPT, Pi);
end