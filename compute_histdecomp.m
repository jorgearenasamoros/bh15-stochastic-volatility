function hd = compute_histdecomp(ma, Ures)
% COMPUTE_HISTDECOMP  Historical decomposition.
%   hd(i,k,t) = contribution of structural shock k to variable i at time t,
%   in the modeled space (here, growth rates):
%       hd(i,k,t) = sum_{j=0}^{t-1} [Psi_j inv(A)]_{ik} * u_{k,t-j}.
%
% Inputs:
%   ma   : n x n x T  unit MA coefficients (ma(:,:,j+1) = Psi_j inv(A)), from structural_ma
%   Ures : T x n      realized structural shocks (row t = u_t' = (A y_t - B x_{t-1})')
[n,~,T] = size(ma);
hd = zeros(n,n,T);
for t = 1:T
    acc = zeros(n,n);
    for j = 0:t-1
        acc = acc + ma(:,:,j+1).*Ures(t-j,:);   % (i,k): Psi_j_invA(i,k) * u_k(t-j)
    end
    hd(:,:,t) = acc;
end
end
