%% COMPUTEBETA_PPT  Computes the density matrix and the optimal value of the SDP involving the PPT relaxation
%   This function has one required input argument:
%       PI: input matrix
%
%   [RHO_PPT, BETA_PPT] = computebeta_PPT(PI) computes the density matrix
%   RHO_PPT and the optimal value of the SDP BETA_PPT
%
%   URL: https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/computebeta_PPT.html
%
%   requires: cvx (http://cvxr.com/cvx), PartialTranspose
%   (http://qetlab.com/PartialTranspose), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function [rho_PPT, beta_PPT] = computebeta_PPT(Pi)
    cvx_begin sdp quiet
        variable rho_PPT(size(Pi)) complex semidefinite;
        maximize HSIP(rho_PPT, Pi)
        subject to
            % density matrix constraint
            trace(rho_PPT) == 1;
            % PPT criterion
            PartialTranspose(rho_PPT, 2) >= 0;
    cvx_end
    beta_PPT = HSIP(rho_PPT, Pi);
end