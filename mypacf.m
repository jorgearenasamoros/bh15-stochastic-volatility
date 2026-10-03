function p = mypacf(x, K)
% MYPACF  Sample partial autocorrelation phi_kk, k=1..K, via Levinson-Durbin.
r   = myacf(x, K);            % rho_1..rho_K
p   = zeros(K,1);
phi = zeros(K,K);
phi(1,1) = r(1);
p(1) = r(1);
v = 1 - r(1)^2;
for k = 2:K
    num = r(k) - sum(phi(k-1,1:k-1)'.*r(k-1:-1:1));
    phi(k,k) = num / max(v,1e-12);
    for j = 1:k-1
        phi(k,j) = phi(k-1,j) - phi(k,k)*phi(k-1,k-j);
    end
    v = v*(1 - phi(k,k)^2);
    p(k) = phi(k,k);
end
end
