% SET_PRIOR  Prior of the BH15 labor-market model (script).
% From setPrior1a.m in the Baumeister and Hamilton (2015) replication files.

% prior for A: Student-t on the demand slope beta (< 0) and the supply slope alpha (> 0)
nA = 2;
cA = [-0.6;0.6];                 % location
sigA = 0.6*ones(nA,1);           % scale
nuA = 3*ones(nA,1);              % degrees of freedom
signA = [-1;1];                  % sign restrictions
longA = [0; 1];                  % equations with a long-run restriction
Ri = [kron(ones(1,nlags),[1 0]) 0];   % demand shock has no long-run effect on employment
Vi = 0.1;
anames = {' \beta'; ' \alpha'};

% prior for D
kappa = 2.0*ones(n,1);

% prior for B (Minnesota-type)
eta = [eye(n) zeros(n,k-n)];
lambda0 = 0.2;
lambda1 = 1.0;
lambda3 = 100.0;
