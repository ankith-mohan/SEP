%%  COMPUTEBETA_R    Compute the density matrix and the optimal value for the SDP involving the Realignment relaxation
%   This function has one required input argument:
%       PI: input matrix
%
%   [RHO_R, BETA_R] = computebeta_r(PI) computes the density matrix RHO_R
%   and the optimal value of the SDP BETA_R
%
%   URL: https://ankith-mohan.github.io/SEP/SDPs/UpperBounds/computebeta_r.html
%
%   requires: cvx (http://cvxr.com/cvx), TraceNorm (http://www.qetlab.com/TraceNorm), Realignment (http://www.qetlab.com/Realignment), HSIP.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


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