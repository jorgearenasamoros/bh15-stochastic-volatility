function [phi, sigma2, acc] = sample_psi_sv(h, phi, a_sig, b_sig, mu_phi, V_phi)
% SAMPLE_PSI_SV  Draw the SV hyperparameters (phi, sigma2) of a log-volatility
%                path, conditional on the path h.
%
%   h_t = phi*h_{t-1} + sigma*eta_t,   eta_t ~ N(0,1),   |phi|<1,
%   h_1 ~ N(0, sigma2/(1-phi^2))   (stationary, zero-mean).
%
% Priors: sigma2 ~ IG(a_sig,b_sig),  phi ~ N(mu_phi,V_phi) truncated to (-1,1).
%
% sigma2 | h,phi is inverse-gamma (conjugate):
%   IG( a_sig + T/2 , b_sig + 0.5*[ sum_{t=2}^T (h_t-phi h_{t-1})^2 + (1-phi^2) h_1^2 ] ).
%
% phi | h,sigma2 is NOT conjugate, because the stationary density of h_1 contributes
%   g(phi) = (1-phi^2)^{1/2} * exp{ -(1-phi^2) h_1^2 / (2 sigma2) }.
% It is drawn by independence Metropolis-Hastings (Kim, Shephard and Chib, 1998). The
% proposal is the Normal implied by the prior and the terms t >= 2, truncated to (-1,1),
%   Vhat  = ( 1/V_phi + sigma2^{-1} sum_{t>=2} h_{t-1}^2 )^{-1}
%   muhat = Vhat*( mu_phi/V_phi + sigma2^{-1} sum_{t>=2} h_{t-1} h_t ),
% drawn exactly by inverting the truncated CDF, and accepted with probability
% min{1, g(phi_new)/g(phi)}. acc = 1 if the proposal was accepted.

h    = h(:);
T    = numel(h);
hlag = h(1:end-1);     % h_{t-1}, t = 2..T
hcur = h(2:end);       % h_t    , t = 2..T

% ---- sigma2 | h, phi   (inverse-gamma) ----
ssr   = sum((hcur - phi*hlag).^2) + (1 - phi^2)*h(1)^2;
shape = a_sig + T/2;
scale = b_sig + 0.5*ssr;
sigma2 = 1 / gamrnd(shape, 1/scale);          % IG(shape,scale)

% ---- phi | h, sigma2   (independence MH with a truncated-Normal proposal) ----
Vhat  = 1 / ( 1/V_phi + sum(hlag.^2)/sigma2 );
muhat = Vhat * ( mu_phi/V_phi + sum(hlag.*hcur)/sigma2 );
phi_prop = rtnorm(muhat, sqrt(Vhat), -1, 1);

logg = @(f) 0.5*log(1 - f^2) - (1 - f^2)*h(1)^2/(2*sigma2);
acc = log(rand) < logg(phi_prop) - logg(phi);
if acc
    phi = phi_prop;
end
end

function x = rtnorm(mu, s, lo, hi)
% Exact draw from N(mu, s^2) truncated to (lo, hi), by inverting the CDF. Uses the
% upper-tail form when the interval lies above the mean, for numerical accuracy.
a = (lo - mu)/s; b = (hi - mu)/s;
if a > 0
    Qa = 0.5*erfc(a/sqrt(2)); Qb = 0.5*erfc(b/sqrt(2));       % upper-tail probabilities
    u  = Qb + rand*(Qa - Qb);
    z  = sqrt(2)*erfcinv(2*u);
else
    Pa = 0.5*erfc(-a/sqrt(2)); Pb = 0.5*erfc(-b/sqrt(2));
    u  = Pa + rand*(Pb - Pa);
    z  = -sqrt(2)*erfcinv(2*u);
end
x = mu + s*z;
x = min(max(x, lo + eps), hi - eps);
end
