function [phi, sigma2] = sample_psi_sv(h, phi, a_sig, b_sig, mu_phi, V_phi)
% SAMPLE_PSI_SV  Draw the SV hyperparameters (phi, sigma2) of the single
%                common log-volatility factor, conditional on the path h.
%
%   h_t = phi*h_{t-1} + sigma*eta_t,   eta_t ~ N(0,1),   |phi|<1,
%   h_1 ~ N(0, sigma2/(1-phi^2))   (stationary, zero-mean).
%
% Conjugate full conditionals:
%   sigma2 | h,phi ~ IG( a_sig + T/2 ,
%                        b_sig + 0.5*[ sum_{t=2}^T (h_t-phi h_{t-1})^2 + (1-phi^2) h_1^2 ] )
%   phi | h,sigma2 ~ N(muhat, Vhat) truncated to (-1,1), with
%       Vhat  = ( 1/V_phi + sigma2^{-1} sum_{t>=2} h_{t-1}^2 )^{-1}
%       muhat = Vhat*( mu_phi/V_phi + sigma2^{-1} sum_{t>=2} h_{t-1} h_t )
%
% Priors: sigma2 ~ IG(a_sig,b_sig),  phi ~ N(mu_phi,V_phi) on (-1,1).

h    = h(:);
T    = numel(h);
hlag = h(1:end-1);     % h_{t-1}, t = 2..T
hcur = h(2:end);       % h_t    , t = 2..T

% ---- sigma2 | h, phi   (inverse-gamma) ----
ssr   = sum((hcur - phi*hlag).^2) + (1 - phi^2)*h(1)^2;
shape = a_sig + T/2;
scale = b_sig + 0.5*ssr;
sigma2 = 1 / gamrnd(shape, 1/scale);          % IG(shape,scale)

% ---- phi | h, sigma2   (Normal truncated to (-1,1)) ----
Vhat  = 1 / ( 1/V_phi + sum(hlag.^2)/sigma2 );
muhat = Vhat * ( mu_phi/V_phi + sum(hlag.*hcur)/sigma2 );

phi_new = muhat + sqrt(Vhat)*randn;
ntry = 0;
while (abs(phi_new) >= 1) && (ntry < 100)      % rejection to enforce stationarity
    phi_new = muhat + sqrt(Vhat)*randn;
    ntry = ntry + 1;
end
if abs(phi_new) < 1
    phi = phi_new;                             % else keep previous phi
end
end
