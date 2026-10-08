function [yall, s, u] = mc_simulate(C, design, rep)
%MC_SIMULATE  Simulates one sample of T = 200 quarters from the BH15 SVAR.
%
%   [yall, s, u] = mc_simulate(C, design, rep)
%
%   A0 y_t = B0 x_{t-1} + u_t,   u_it = sqrt(d0_i * s_it) * z_it,   z_it ~ N(0,1) iid,
%   x_{t-1} = [y_{t-1}' ... y_{t-8}' 1]',
%
%   with the true parameters in C (calib.mat). The design sets the variance scale s_it:
%     'dgp1'  : constant volatility, s_it = 1;
%     'dgp2k' : common factor, s_it = C.s2k(t) for both shocks;
%     'dgp3'  : shock-specific volatility, s_it = C.s3(t,i).
%   In every design mean_t(s_it) = 1, so the average variance of shock i is d0_i.
%   The normals z depend only on the replication (seed C.seed_z + rep), so the three designs
%   use the same random numbers. Initial conditions: the 8 observed quarters before 1970:Q1 (C.y0).
%
%   Output: yall (208 x 2: 8 initial and 200 simulated quarters; the estimation sample is rows
%           9 to 208), s (200 x 2) variance scales, u (200 x 2) structural shocks.
T = C.T; n = 2; p = 8;
rng(C.seed_z + rep, 'twister');
z = randn(T, n);
switch design
    case 'dgp1'
        s = ones(T, n);
    case 'dgp2k'
        s = repmat(C.s2k, 1, n);
    case 'dgp3'
        s = C.s3;
    otherwise
        error('mc_simulate: unknown design %s', design);
end
u   = sqrt(C.d0 .* s) .* z;              % T x n (d0 is 1 x n)
ept = u / C.A0';                         % reduced-form errors (A0 \ u_t')'
y = [C.y0; zeros(T, n)];
for t = 1:T
    x = [reshape(y(t+p-1:-1:t, :)', [], 1); 1];   % [w(-1) n(-1) ... w(-8) n(-8) 1]
    y(t+p, :) = (C.Phi0*x)' + ept(t, :);
end
yall = y;
end
