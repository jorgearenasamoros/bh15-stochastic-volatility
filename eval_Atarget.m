function [ptarget, zeta, mstar, tau] = eval_Atarget(vecA,cA,sigA,nuA,signA, ...
        T,omegahat,Sstar,ytilde,yxtilde,xtildei,kappa,kappastar,longA,Ri,Vi)
% EVAL_ATARGET  Log of the (B,D)-marginalized posterior kernel for A in the
%               BH15 conjugate SVAR, plus the by-products needed to draw (b_i,d_ii).
%
% This is the target used inside the RW-Metropolis loop of the
% original BH15 replication code (main_labor.m), refactored into a function so it can be evaluated
% at both the current and proposed A, and so the (b_i,d_ii) draw can reuse
% zeta/mstar/tau.  The ONLY change for stochastic volatility is that the
% sufficient statistics (ytilde,yxtilde,xtildei) are passed in already
% volatility-standardized (weights e^{-h_t}); the algebra is unchanged.
%
% Inputs:
%   vecA            : free elements of A
%   cA,sigA,nuA,signA : Student-t prior parameters for A (see logP.m)
%   T               : sample size
%   omegahat        : OLS residual covariance (constant; supplies |det A|^T via det term)
%   Sstar           : univariate residual covariance (prior centring of d_ii)
%   ytilde          : YY'*W*YY + eta*Mtildeinv*eta'        (n x n)
%   yxtilde         : YY'*W*XX + eta*Mtildeinv             (n x k)
%   xtildei         : inv(XX'*W*XX + Mtildeinv [+ Ri'Ri/Vi]) per equation = Mstar_i (k x k x n)
%   kappa,kappastar : IG shape for d_ii, prior and posterior
%   longA,Ri,Vi     : long-run prior restriction (per equation)
%
% Outputs:
%   ptarget : log posterior kernel at vecA (up to A-independent constants)
%   zeta    : n x 1, with tau_i + zeta_i/2 = tau*_i (posterior IG scale)
%   mstar   : n x k, row i = bhat_i'  (= yxtildei*Mstar_i)
%   tau     : n x 1 prior IG scale, kappa.*diag(A*Sstar*A')

A = setA(vecA);
n = size(omegahat,1);
k = size(yxtilde,2);

Q   = A*omegahat*A';                 % det term -> (T/2)log det Q = T*log|det A| + const
tau = kappa.*diag(A*Sstar*A');       % prior IG scale (depends on A, NOT on volatility)

zeta  = zeros(n,1);
mstar = zeros(n,k);
i = 0;
while i < n
    i = i + 1;
    ytildei  = A(i,:)*ytilde*A(i,:)';
    yxtildei = A(i,:)*yxtilde;
    if longA(i) == 1
        ri = A(2,1);
        ytildei  = ytildei  + ri^2/Vi;
        yxtildei = yxtildei + ri*Ri/Vi;
    end
    zeta(i)    = ytildei - yxtildei*xtildei(:,:,i)*yxtildei';
    mstar(i,:) = yxtildei*xtildei(:,:,i);
end

ptarget = logP(vecA,cA,sigA,nuA,signA) + (T/2)*log(det(Q)) ...
        - kappastar'*log((2*tau/T) + zeta/T) + kappa'*log(tau);
end
