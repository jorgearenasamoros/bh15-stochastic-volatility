function out = mc_estimate(yall, tstart, tend, SV_on, opt)
%MC_ESTIMATE  Gibbs sampler of the BH15 labor-market SVAR, with or without the common
%   stochastic-volatility factor, written as a function for the Monte Carlo.
%
%   out = mc_estimate(yall, tstart, tend, SV_on, opt)
%
%   Same model, priors and blocks as gibbs_sv.m:
%     1. A | Y,H           random-walk Metropolis on theta = log(signA.*a), with the Jacobian, on
%                          the (B,D)-marginal posterior; the proposal covariance is adapted
%                          (Haario et al., 2001) during the burn-in only, every 500 iterations from
%                          iteration 2,000;
%     2. (b_i,d_ii) | A,H  conjugate Normal-inverse-gamma draws;
%     3. h | A,B,D         Kim-Shephard-Chib mixture and scalar FFBS (sample_h_common.m), with h
%                          recentred to mean zero at every sweep;
%     4. (phi,sigma2) | h  sample_psi_sv.m.
%   Differences from gibbs_sv.m, none of which changes the target distribution:
%     - the chain starts at the posterior mode of A (fminunc on post_val.m, the log posterior with
%       B and D integrated out) if it satisfies the sign restrictions, otherwise at the prior
%       location cA;
%     - the seed is set with rng(seed,'twister');
%     - with SV_on = 0 the log target at the current A is not recomputed at the start of each
%       sweep (with h = 0 the sufficient statistics do not change, and no random numbers are used);
%     - B, 1./d and h are kept every opt.thin draws; A, phi and sigma2 at every draw.
%
%   Inputs:
%     yall         : Tall x 2 [wage growth, employment growth] in 100*dlog, with at least 8
%                    quarters before tstart
%     tstart, tend : rows of yall in the estimation sample
%     SV_on        : 1 = common stochastic-volatility factor, 0 = homoskedastic
%     opt          : struct with NIT (total iterations), NB (burn-in), seed, thin (default 5)
%
%   Output (struct): a (2 x NIT-NB draws of [beta; alpha]), B (2 x 17 x nthin), invd (2 x nthin),
%     h (T x nthin), phi and sig2 (1 x NIT-NB) with SV, acc_all and acc_post (acceptance rates over
%     all iterations and after the burn-in), a_start (starting point), secs.
%
%   Requires the Optimization Toolbox (fminunc) and the Statistics and Machine Learning Toolbox.

if ~isfield(opt,'thin') || isempty(opt.thin), opt.thin = 5; end
t0 = tic;
ndraws = opt.NIT; nburn = opt.NB;

%% ---------------- data ----------------
nlags = 8;
YY = yall(tstart:tend,:);
[T, n] = size(YY);
XX = zeros(T, n*nlags + 1);
for l = 1:nlags
    XX(:,(l-1)*n+(1:n)) = yall(tstart-l:tend-l,:);    % [w(-1) n(-1) w(-2) n(-2) ... 1]
end
XX(:,end) = 1;
k = size(XX,2);
omegahat = (YY'*YY - YY'*XX*inv(XX'*XX)*XX'*YY)/T;    %#ok<MINV>

%% ---------------- priors (as in set_prior.m and gibbs_sv.m) ----------------
nA = 2;
cA = [-0.6; 0.6];  sigA = 0.6*ones(nA,1);  nuA = 3*ones(nA,1);  signA = [-1; 1];
longA = [0; 1];
Ri = [kron(ones(1,nlags),[1 0]) 0];  Vi = 0.1;
kappa = 2.0*ones(n,1);
eta = [eye(n) zeros(n,k-n)];
lambda0 = 0.2;  lambda1 = 1.0;  lambda3 = 100.0;
% stochastic volatility
offset = 1e-6;                      % avoids log(0) in log(eps^2)
a_sig = 5;  b_sig = 0.2;            % sigma2 ~ IG(5, 0.2)
mu_phi = 0.95;  V_phi = 0.10^2;     % phi ~ N(0.95, 0.1^2) truncated to (-1,1)
phi = 0.90;  sigma2 = 0.05;
h = zeros(T,1);
[qK,mK,v2K] = ksc7_mixture();

%% ---------------- Sstar: residual covariance of univariate AR(8) regressions ----------------
e = zeros(T,n);
for i = 1:n
    ylags = zeros(T, nlags+1);
    for l = 1:nlags, ylags(:,l) = yall(tstart-l:tend-l,i); end
    ylags(:,end) = 1;
    e(:,i) = yall(tstart:tend,i) - ylags*inv(ylags'*ylags)*ylags'*yall(tstart:tend,i); %#ok<MINV>
end
Sstar = e'*e/T;

%% ---------------- Minnesota prior ----------------
mv1 = (1:nlags)'.^(-2*lambda1);
mv2 = 1./diag(Sstar);
mv3 = kron(mv1,mv2);
mv3 = lambda0^2*[mv3; lambda3^2];
mv3 = 1./mv3;
Mtildeinv = diag(mv3);
kappastar = kappa + (T/2);
priorYY = eta*Mtildeinv*eta';
priorYX = eta*Mtildeinv;

% homoskedastic sufficient statistics
ytilde0  = YY'*YY + priorYY;
yxtilde0 = YY'*XX + priorYX;
xtildei0 = zeros(k,k,n);
for i = 1:n
    M = XX'*XX + Mtildeinv;
    if longA(i) == 1, M = M + Ri'*Ri/Vi; end
    xtildei0(:,:,i) = inv(M);
end

rng(opt.seed, 'twister');

%% ---------------- starting point: posterior mode of A ----------------
A_params = [cA sigA nuA signA];
f_anon = @(theta) post_val(theta,A_params,longA,kappa,T,omegahat,Sstar,ytilde0,xtildei0,yxtilde0,Ri,Vi);
theta_max = fminunc(f_anon, cA, optimset('Display','off'));
if min(sign(theta_max).*signA) >= 0
    a_old = theta_max;
else
    a_old = cA;
end
a_start = a_old;

%% ---------------- storage ----------------
nsave  = ndraws - nburn;
nthin  = numel(opt.thin:opt.thin:nsave);
a_post = zeros(nA, nsave);
B_post = zeros(n, k, nthin);
invd_post = zeros(n, nthin);
if SV_on
    h_post = zeros(T, nthin);  phi_post = zeros(1, nsave);  sig2_post = zeros(1, nsave);
end

%% ---------------- random walk on log|a| for the sign-restricted elements ----------------
sc  = signA + (signA==0);
isc = signA ~= 0;
a2t  = @(a) [log(sc(isc).*a(isc)); a(~isc)];          % a -> theta
logJ = @(a) sum(log(abs(a(isc))));                     % log |da/dtheta|
Lprop = 0.3*eye(nA);                                   % initial scale, adapted in the burn-in
th_hist = zeros(nA, max(nburn,1));
naccept = 0;  nacc_post = 0;

ytilde = ytilde0;  yxtilde = yxtilde0;  xtildei = xtildei0;
jt = 0;
[pt_old, zeta_old, mstar_old, tau_old] = eval_Atarget(a_old,cA,sigA,nuA,signA, ...
        T,omegahat,Sstar,ytilde,yxtilde,xtildei,kappa,kappastar,longA,Ri,Vi);

for count = 1:ndraws
    % ---- volatility-standardized sufficient statistics given the current path h ----
    if SV_on
        wv = exp(-h);
        ytilde  = (YY.*wv)'*YY + priorYY;
        yxtilde = (YY.*wv)'*XX + priorYX;
        XXw = (XX.*wv)'*XX;
        for i = 1:n
            M = XXw + Mtildeinv;
            if longA(i) == 1, M = M + Ri'*Ri/Vi; end
            xtildei(:,:,i) = inv(M);
        end
        [pt_old, zeta_old, mstar_old, tau_old] = eval_Atarget(a_old,cA,sigA,nuA,signA, ...
                T,omegahat,Sstar,ytilde,yxtilde,xtildei,kappa,kappastar,longA,Ri,Vi);
    end

    % ---- Step 1: A | Y,H ----
    th_old = a2t(a_old);
    a_new  = local_t2a(th_old + Lprop*randn(nA,1), sc, isc);
    [pt_new, zeta_new, mstar_new, tau_new] = eval_Atarget(a_new,cA,sigA,nuA,signA, ...
            T,omegahat,Sstar,ytilde,yxtilde,xtildei,kappa,kappastar,longA,Ri,Vi);
    if log(rand) <= (pt_new + logJ(a_new)) - (pt_old + logJ(a_old))
        a_old = a_new;  pt_old = pt_new;  zeta_old = zeta_new;  mstar_old = mstar_new;  tau_old = tau_new;
        naccept = naccept + 1;
        if count > nburn, nacc_post = nacc_post + 1; end
    end
    if count <= nburn
        th_hist(:,count) = a2t(a_old);
        if count >= 2000 && mod(count,500) == 0
            Cv = cov(th_hist(:,ceil(count/2):count)');
            Lprop = chol((2.38^2/nA)*Cv + 1e-8*eye(nA))';
        end
    end

    % ---- Step 2: (b_i,d_ii) | A,H ----
    A    = setA(a_old);
    invd = zeros(n,1);
    Bd   = zeros(n,k);
    for i = 1:n
        taustar_i = tau_old(i) + zeta_old(i)/2;
        invd(i)   = gamrnd(kappastar(i), 1/taustar_i);
        d_ii      = 1/invd(i);
        bhat      = (mstar_old(i,:))';
        Bd(i,:)   = (bhat + sqrt(d_ii)*(chol(xtildei(:,:,i))'*randn(k,1)))';
    end

    % ---- Steps 3-4: volatility path and its hyperparameters ----
    if SV_on
        Uresid  = YY*A' - XX*Bd';
        EPS     = Uresid ./ sqrt((1./invd)');
        eps2log = log(EPS.^2 + offset);
        h       = sample_h_common(eps2log,h,phi,sigma2,qK,mK,v2K);
        h       = h - mean(h);          % the levels of h_t and d_ii are not separately identified
        [phi,sigma2] = sample_psi_sv(h,phi,a_sig,b_sig,mu_phi,V_phi);
    end

    % ---- store ----
    if count > nburn
        j = count - nburn;
        a_post(:,j) = a_old;
        if SV_on, phi_post(j) = phi; sig2_post(j) = sigma2; end
        if mod(j, opt.thin) == 0
            jt = jt + 1;
            B_post(:,:,jt) = Bd;  invd_post(:,jt) = invd;
            if SV_on, h_post(:,jt) = h; end
        end
    end
end

out.a = a_post;  out.B = B_post;  out.invd = invd_post;
if SV_on
    out.h = h_post;  out.phi = phi_post;  out.sig2 = sig2_post;
end
out.acc_all  = naccept/ndraws;
out.acc_post = nacc_post/max(nsave,1);
out.a_start  = a_start;
out.thin = opt.thin;  out.NIT = ndraws;  out.NB = nburn;  out.seed = opt.seed;  out.SV_on = SV_on;
out.secs = toc(t0);
end

function a = local_t2a(t, sc, isc)
% theta -> a: exponential for the sign-restricted elements
a = zeros(numel(sc),1);
a(isc)  = sc(isc).*exp(t(1:nnz(isc)));
a(~isc) = t(nnz(isc)+1:end);
end

function [q, m, v2] = ksc7_mixture()
% Kim, Shephard and Chib (1998), Table 4: 7-component normal mixture for log(chi^2_1), the same
% constants as ksc7.m. The component means are the Table 4 values plus E[log chi^2_1] =
% psi(1/2) + log(2), evaluated here in double precision (ksc7.m uses the rounded -1.2703628; the
% difference, 4.5e-8, is immaterial, and this form reproduces the stored Monte Carlo estimates
% bit for bit).
q     = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];
m_tab = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];
v2    = [ 5.79596;  2.61369;  5.17950; 0.16735; 0.64009; 0.34023; 1.26261];
m     = m_tab + (psi(0.5) + log(2));
end
