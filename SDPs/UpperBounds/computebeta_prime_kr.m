%%  COMPUTEBETA_PRIME_KR    Computes the density matrix and the optimal value for the SDP involving the Bosonic Extensions as well as the Realignment relaxations
%   This function has four required arguments:
%       PI: input matrix
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       K: integer that describes the number of bosonic extensions
%       PPT: 1 to impose PPT criterion
%            0 otherwise
%
%   [RHO_PRIME_KR, BETA_PRIME_KR] = computebeta_prime_kr(PI, DIM_B, K, PPT)
%   computes the density matrix RHO_PRIME_KR and the optiaml value of the
%   SDP BETA_PRIME_KR
%
%   URL:
%   https://ankith-mohan.github.io/SEP/helpers/computebeta_prime_kr.html
%
%   requires: cvx (http://cvxr.com/cvx), SymmetricExtension
%   (http://www.qetlab.com/SymmetricExtension), TraceNorm
%   (http://www.qetlab.com/TraceNorm), Realignment
%   (http://www.qetlab.com/Realignment), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr(Pi, dim_B, k, ...
                                                            PPT)
    cvx_begin sdp quiet
        variable rho_prime_kr(size(Pi)) complex semidefinite;
        maximize HSIP(rho_prime_kr, Pi)
        subject to
            trace(rho_prime_kr) == 1; % density matrix
            SymmetricExtension(rho_prime_kr, k, dim_B, PPT, 1, eps^(1/4)); 
            % Bosonic symmetric extension constraint
            TraceNorm(Realignment(rho_prime_kr)) <= 1; 
            % realignment criterion
    cvx_end
    beta_prime_kr = HSIP(rho_prime_kr, Pi);
end