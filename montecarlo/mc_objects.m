function [Z, names] = mc_objects(a, B, d)
%MC_OBJECTS  The 22 objects evaluated in the Monte Carlo, for each posterior draw.
%
%   [Z, names] = mc_objects(a, B, d)
%     a : 2 x S draws of [beta; alpha];  B : 2 x 17 x S;  d : 2 x S structural variances d_ii
%     Z : S x 22, columns in the order of names
%
%   Columns:
%     1-2    beta, alpha (free elements of A = [-beta 1; -alpha 1])
%     3-6    sums of the 8 lag coefficients of the reduced form Phi = A\B: wage equation (own,
%            cross), employment equation (own, cross)
%     7-8    constants of the reduced form (wage, employment)
%     9-20   responses to a unit structural shock, cumulated to log levels (x100), at 0, 4 and 20
%            quarters: wage <- demand, wage <- supply, employment <- demand, employment <- supply
%     21-22  forecast-error variance shares of the growth rates at 16 quarters: demand in wages,
%            supply in employment
%   The impact responses to a supply shock (columns 12 and 18) are implied by those to a demand
%   shock (wage <- supply = -(wage <- demand), employment <- supply = 1 - (employment <- demand)),
%   so the summaries over objects use the other 20 columns.
nlags = 8; n = 2; hmax = 21; Hf = 16; hix = [1 5 21];
S = size(a,2);
Z = zeros(S, 22);
iw = 1:2:2*nlags; in = 2:2:2*nlags;
for s = 1:S
    A = setA(a(:,s)); Bs = B(:,:,s);
    Phi = A\Bs;
    ma  = structural_ma(A, Bs, nlags, n, hmax);
    cma = cumsum(ma, 3);
    c = zeros(n,n);
    for hh = 1:Hf, c = c + (ma(:,:,hh).^2).*(d(:,s)'); end
    fe = c./sum(c,2);
    Z(s,:) = [a(1,s) a(2,s) ...
        sum(Phi(1,iw)) sum(Phi(1,in)) sum(Phi(2,in)) sum(Phi(2,iw)) Phi(1,end) Phi(2,end) ...
        squeeze(cma(1,1,hix))' squeeze(cma(1,2,hix))' squeeze(cma(2,1,hix))' squeeze(cma(2,2,hix))' ...
        fe(1,1) fe(2,2)];
end
names = {'beta','alpha', ...
    'LR_w_own','LR_w_cross','LR_n_own','LR_n_cross','const_w','const_n', ...
    'irf_w_dem_h0','irf_w_dem_h4','irf_w_dem_h20', ...
    'irf_w_sup_h0','irf_w_sup_h4','irf_w_sup_h20', ...
    'irf_n_dem_h0','irf_n_dem_h4','irf_n_dem_h20', ...
    'irf_n_sup_h0','irf_n_sup_h4','irf_n_sup_h20', ...
    'fevd_w_dem_h16','fevd_n_sup_h16'};
end
