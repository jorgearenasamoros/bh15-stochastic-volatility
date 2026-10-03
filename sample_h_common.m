function h = sample_h_common(eps2log, h, phi, sigma2, q, m, v2)
% SAMPLE_H_COMMON  Draw the single common log-volatility path h_{1:T}.
%
% Model (single common factor):
%   measurement (n obs per period, all loading on the scalar state):
%       log(eps_it^2) = h_t + m_{s_it} + sqrt(v2_{s_it})*zeta_it,   i = 1..n
%   state:
%       h_t = phi*h_{t-1} + sigma*eta_t,   h_1 ~ N(0, sigma2/(1-phi^2)).
%
% Steps:
%   (1) sample KSC mixture indicators s_it (vectorized over the T x n grid);
%   (2) collapse the n conditionally-Gaussian measurements at each t into a
%       single equivalent observation (precision-weighted) -> scalar state space;
%   (3) scalar forward-filter / backward-sample (FFBS).
%
% Inputs:
%   eps2log : T x n matrix of log(eps_it^2 + offset), eps_it = u_it/sqrt(d_ii)
%   h       : current T x 1 path (used only to draw the indicators)
%   phi,sigma2 : AR(1) parameters
%   q,m,v2  : KSC 7-component constants (each 7x1), from ksc7()
% Output:
%   h       : new T x 1 path.

[T, n] = size(eps2log);
K = numel(q);

% ---- (1) mixture indicators s_it (vectorized) ----
r = eps2log - h;                       % T x n residual about current state (implicit expansion)
logw = zeros(T, n, K);
for j = 1:K
    logw(:,:,j) = log(q(j)) - 0.5*log(v2(j)) - 0.5*((r - m(j)).^2)/v2(j);
end
logw = logw - max(logw, [], 3);        % stabilize
w    = exp(logw);
w    = w ./ sum(w, 3);
cw   = cumsum(w, 3);
U    = rand(T, n);
S    = sum(cw < U, 3) + 1;             % T x n indicator in 1..K
S    = min(S, K);

ms = m(S);                             % T x n component means
vs = v2(S);                            % T x n component variances

% ---- (2) collapse n Gaussian measurements per period into one ----
obs  = eps2log - ms;                   % T x n  (each ~ N(h_t, vs))
Rinv = sum(1./vs, 2);                  % T x 1  combined precision
ybar = sum(obs./vs, 2) ./ Rinv;        % T x 1  combined observation
R    = 1 ./ Rinv;                      % T x 1  combined measurement variance

% ---- (3a) scalar Kalman forward filter ----
a_f = zeros(T,1);                      % filtered mean E[h_t | y_{1:t}]
P_f = zeros(T,1);                      % filtered variance
a_pred = 0;                            % prior mean of h_1
P_pred = sigma2/(1 - phi^2);           % stationary prior variance of h_1
for t = 1:T
    if t > 1
        a_pred = phi*a_f(t-1);
        P_pred = phi^2*P_f(t-1) + sigma2;
    end
    Kt     = P_pred / (P_pred + R(t));
    a_f(t) = a_pred + Kt*(ybar(t) - a_pred);
    P_f(t) = (1 - Kt)*P_pred;
end

% ---- (3b) backward sampling ----
h    = zeros(T,1);
h(T) = a_f(T) + sqrt(P_f(T))*randn;
for t = T-1:-1:1
    Ppred = phi^2*P_f(t) + sigma2;     % P_{t+1|t}
    J     = phi*P_f(t)/Ppred;
    mu    = a_f(t) + J*(h(t+1) - phi*a_f(t));
    V     = P_f(t) - J*phi*P_f(t);     % = P_f(t) - (phi P_f(t))^2/Ppred
    h(t)  = mu + sqrt(max(V,0))*randn;
end
end
