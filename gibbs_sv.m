% GIBBS_SV  Gibbs sampler for the BH15 labor-market SVAR with a single common
% stochastic-volatility factor (script, called by run_chain.m).
%
% Extends main_labor.m of the Baumeister and Hamilton (2015) replication files:
%     A y_t = B x_{t-1} + u_t,   u_t ~ N(0, e^{h_t} D),
%     h_t = phi h_{t-1} + sigma eta_t,   h_1 ~ N(0, sigma2/(1-phi^2)).
%
% Gibbs sweep:
%   1. A | Y,H          random-walk Metropolis on the (B,D)-marginal posterior,
%                       computed with volatility-standardized data. The walk runs
%                       on theta = log(signA.*a) for sign-restricted elements, with
%                       the Jacobian, so every proposal satisfies the restrictions.
%                       The proposal covariance is adapted during the burn-in only.
%   2. (b_i,d_ii) | A,H conjugate Normal-inverse-gamma draws
%   3. h_{1:T} | A,B,D  Kim-Shephard-Chib mixture and scalar FFBS (sample_h_common)
%   4. (phi,sigma2) | H inverse-gamma draw for sigma2, Metropolis-Hastings step for phi (sample_psi_sv)
% With SV_on = 0, h_t = 0 and the sampler is the homoskedastic BH15 model.
%
% Set by the caller: SV_on, ndraws, nburn, SEED, A_START.

%% data
read_data;                                   % YY, yall, tstart, tend, nlags, time
T = size(YY,1);
n = size(YY,2);
XX = yall(tstart-1:tend-1,:);
for ilags = 2:nlags
    XX = [XX yall(tstart-ilags:tend-ilags,:)]; %#ok<AGROW>
end
XX = [XX ones(T,1)];
k = size(XX,2);
omegahat = (YY'*YY - YY'*XX*inv(XX'*XX)*XX'*YY)/T;  %#ok<MINV>

%% priors
set_prior;

% SV priors and initial values
offset = 1e-6;                    % avoids log(0) in log(eps^2)
a_sig  = 5;    b_sig = 0.2;       % sigma2 ~ IG(a_sig,b_sig), prior mean 0.05
mu_phi = 0.95; V_phi = 0.10^2;    % phi ~ N(mu_phi,V_phi) truncated to (-1,1)
phi    = 0.90; sigma2 = 0.05;     % initial hyperparameters
h      = zeros(T,1);              % initial path
[qK,mK,v2K] = ksc7();

% univariate AR residual covariance Sstar (centres the prior of d_ii)
e = zeros(T,n);
for i = 1:n
    ylags = yall(tstart-1:tend-1,i);
    for ilags = 2:nlags
        ylags = [ylags yall(tstart-ilags:tend-ilags,i)]; %#ok<AGROW>
    end
    ylags = [ylags ones(T,1)]; %#ok<AGROW>
    e(:,i) = yall(tstart:tend,i) - ylags*inv(ylags'*ylags)*ylags'*yall(tstart:tend,i); %#ok<MINV>
end
Sstar = e'*e/T;

% Minnesota prior precision
mv1 = (1:nlags)'.^(-2*lambda1);
mv2 = 1./diag(Sstar);
mv3 = kron(mv1,mv2);
mv3 = lambda0^2*[mv3; lambda3^2];
Mtildeinv = diag(1./mv3);

kappastar = kappa + (T/2);
priorYY = eta*Mtildeinv*eta';
priorYX = eta*Mtildeinv;

% homoskedastic sufficient statistics (used when SV_on = 0)
ytilde  = YY'*YY + priorYY;
yxtilde = YY'*XX + priorYX;
xtildei = zeros(k,k,n);
for i = 1:n
    M = XX'*XX + Mtildeinv;
    if longA(i) == 1, M = M + Ri'*Ri/Vi; end
    xtildei(:,:,i) = inv(M);
end

rand('seed',SEED);  %#ok<RAND>
randn('seed',SEED); %#ok<RAND>

%% storage
nsave     = ndraws - nburn;
a_post    = zeros(nA,nsave);
invd_post = zeros(n,nsave);      % 1./d_ii
B_post    = zeros(n,k,nsave);
h_post    = zeros(T,nsave);
phi_post  = zeros(1,nsave);
sig2_post = zeros(1,nsave);

%% random walk on log scale for the sign-restricted elements of A
a_old = A_START(:);
isc  = signA ~= 0;
sc   = signA + (signA==0);
a2t  = @(a) isc.*log(abs(a) + (~isc)) + (~isc).*a;       % a -> theta
t2a  = @(t) isc.*sc.*exp(t.*isc) + (~isc).*t;            % theta -> a
logJ = @(a) sum(log(abs(a(isc))));                       % log |da/dtheta|
Lprop = 0.3*eye(nA);                                     % initial scale, adapted in the burn-in
th_hist = zeros(nA, nburn);
naccept = 0;

if SV_on, modeltag = 'SV'; else, modeltag = 'homoskedastic'; end
fprintf('Running %s sampler: %d draws (%d burn-in), seed %d\n', modeltag, ndraws, nburn, SEED);

for count = 1:ndraws

    % volatility-standardized sufficient statistics given the current path h
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
    end
    [pt_old, zeta_old, mstar_old, tau_old] = eval_Atarget(a_old,cA,sigA,nuA,signA, ...
            T,omegahat,Sstar,ytilde,yxtilde,xtildei,kappa,kappastar,longA,Ri,Vi);

    % Step 1: A | Y,H
    a_new = t2a(a2t(a_old) + Lprop*randn(nA,1));
    [pt_new, zeta_new, mstar_new, tau_new] = eval_Atarget(a_new,cA,sigA,nuA,signA, ...
            T,omegahat,Sstar,ytilde,yxtilde,xtildei,kappa,kappastar,longA,Ri,Vi);
    if log(rand) <= (pt_new + logJ(a_new)) - (pt_old + logJ(a_old))
        a_old = a_new;  zeta_old = zeta_new;  mstar_old = mstar_new;  tau_old = tau_new;
        naccept = naccept + 1;
    end
    if count <= nburn                                    % adaptation (Haario et al., 2001)
        th_hist(:,count) = a2t(a_old);
        if count >= 2000 && mod(count,500) == 0
            C = cov(th_hist(:,ceil(count/2):count)');
            Lprop = chol((2.38^2/nA)*C + 1e-8*eye(nA))';
        end
    end

    % Step 2: (b_i,d_ii) | A,H
    A    = setA(a_old);
    invd = zeros(n,1);
    Bd   = zeros(n,k);
    for i = 1:n
        taustar_i = tau_old(i) + zeta_old(i)/2;
        invd(i)   = gamrnd(kappastar(i), 1/taustar_i);
        bhat      = (mstar_old(i,:))';
        Bd(i,:)   = (bhat + sqrt(1/invd(i))*(chol(xtildei(:,:,i))'*randn(k,1)))';
    end

    % Steps 3-4: volatility path and its hyperparameters
    if SV_on
        EPS     = (YY*A' - XX*Bd') ./ sqrt((1./invd)');
        h       = sample_h_common(log(EPS.^2 + offset),h,phi,sigma2,qK,mK,v2K);
        h       = h - mean(h);   % the levels of h_t and d_ii are not separately identified
        [phi,sigma2] = sample_psi_sv(h,phi,a_sig,b_sig,mu_phi,V_phi);
    end

    if count > nburn
        j = count - nburn;
        a_post(:,j)    = a_old;
        invd_post(:,j) = invd;
        B_post(:,:,j)  = Bd;
        h_post(:,j)    = h;
        phi_post(j)    = phi;
        sig2_post(j)   = sigma2;
    end
    if mod(count,10000) == 0
        fprintf('  draw %6d / %d   (acceptance %.3f)\n', count, ndraws, naccept/count);
    end
end
acceptance_ratio = naccept/ndraws;
