function ma = structural_ma(A,B,nlags,n,H)
% STRUCTURAL_MA  Structural moving-average coefficients of the SVAR.
%   ma(:,:,h) = Psi_{h-1} * inv(A) = response of y to a UNIT structural shock,
%   h = 1 (impact) ... H. Psi_j are the reduced-form MA coefficients implied by
%   the reduced-form lag matrices Phi = inv(A)*B.
%
% Inputs: A (n x n), B (n x k, k = n*nlags+1), nlags, n, H (max horizon).
invA = inv(A);
Phi  = invA*B;                       % n x k reduced-form coefficients
psi  = zeros(n,n,H);
psi(:,:,1) = eye(n);                 % Psi_0
for h = 2:H
    S = zeros(n,n);
    for l = 1:min(h-1,nlags)
        S = S + Phi(:,(l-1)*n+1:l*n)*psi(:,:,h-l);
    end
    psi(:,:,h) = S;
end
ma = zeros(n,n,H);
for h = 1:H
    ma(:,:,h) = psi(:,:,h)*invA;     % unit structural IRF
end
end
