function [q,m,v2] = ksc7()
% KSC7  Kim, Shephard & Chib (1998) 7-component Gaussian mixture that
%       approximates the distribution of log(chi-square_1).
%
% Used to linearize the stochastic-volatility measurement equation
%       log(eps_it^2) = h_t + log(xi_it^2),   xi_it ~ N(0,1),
% where log(xi_it^2) ~ log(chi^2_1) is replaced by a mixture of 7 normals
% with indicator s_it in {1,...,7}:
%       log(eps_it^2) = h_t + m_{s_it} + sqrt(v2_{s_it})*zeta_it.
%
% Outputs (each 7x1):
%   q  : component probabilities (sum to 1)
%   m  : component means, including E[log chi^2_1] = psi(1/2) + log(2) = -1.2703628
%   v2 : component variances
%
% The constants in Table 4 of KSC (1998) are centred, sum(q.*m_table) = 0; the component
% means of the mixture for log(chi^2_1) are m_table + E[log chi^2_1] = m_table - 1.2703628.
% The offset matters here because h_t is recentred to mean zero at every sweep, so its
% level cannot absorb a shift in the means.
%
% Reference: Kim, Shephard & Chib (1998, RES), Table 4.
q      = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];
m_tab  = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];   % Table 4 (centred)
m      = m_tab - 1.2703628;                     % add E[log chi^2_1]
v2     = [ 5.79596;  2.61369;  5.17950; 0.16735; 0.64009; 0.34023; 1.26261];
end
