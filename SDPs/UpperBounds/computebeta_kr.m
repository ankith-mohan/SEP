%% COMPUTEBETA_KR   Computes the density matrix and the optimal value for the SDP involving both the Symmetric Extensions as well as the Realignment relaxations
%   This function has four required input arguments:
%       PI: input matrix
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       K: integer that describes the number of symmetric extensions
%       PPT: 1 to impose PPT criterion
%            0 otherwise
%
%   [RHO_KR, BETA_KR] = computebeta_kr(PI, DIM_B, K, PPT) computes the
%   density matrix RHO_KR and the optimal value of the SDP BETA_KR
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/computebeta_kr.html
%
%   requires: cvx (http://cvxr.com/cvx/), SymmetricExtension
%   (http://www.qetlab.com/SymmetricExtension), TraceNorm (http://www.qetlab.com/TraceNorm), Realignment (http://www.qetlab.com/Realignment), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_kr, beta_kr] = computebeta_kr(Pi, dim_B, k, PPT)
    cvx_begin sdp quiet
        variable rho_kr(size(Pi)) complex semidefinite;
        maximize HSIP(rho_kr, Pi)
        subject to
            % density matrix constraint
            trace(rho_kr) == 1;
            % Symmetric extension constraint
            SymmetricExtension(rho_kr, k, dim_B, PPT, 0, eps^(1/4));
            % realignment constraint
            TraceNorm(Realignment(rho_kr)) <= 1;
    cvx_end
    beta_kr = HSIP(rho_kr, Pi);
end