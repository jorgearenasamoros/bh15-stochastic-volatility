function fig_loglik(outfile)
% FIG_LOGLIK  Figure 1 of the post: contours of the log likelihood of the SV model for the supply
% slope alpha and the demand slope beta, given the posterior mean of the volatility path (the
% design of the lower panel of Figure 4 in Baumeister and Hamilton, 2015).
%
% With A = [-beta 1; -alpha 1] and B and D concentrated out, the log likelihood is
%     T log|det A| - (T/2) log det[diag(A Omega A')],    det A = alpha - beta.
% Given the volatility path, the standardized data e^{-h_t/2} y_t have a likelihood of this form
% (the Jacobian of the standardization does not depend on A, B or D). Because the factor is common
% to both equations, B is concentrated out by weighted least squares with weights e^{-h_t}, and
%     Omega* = (1/T) sum_t e^{-hbar_t} e_t e_t',
% where e_t are the weighted least squares residuals and hbar_t is the posterior mean of h_t over
% the four SV chains (results/chain_sv_*.mat, every 10th draw). Contour lines are drawn every 10
% points, from 10 to 100 below the maximum. The black curve is the set of maximum likelihood
% estimates, beta(alpha) = (w22 - alpha w12)/(w12 - alpha w11) (BH15, eq. 51), the pairs that make
% A Omega* A' diagonal; the likelihood is the same at every point of it.
%
% fig_loglik(outfile) saves the figure to outfile (default figures/loglik-sv.png).
rootdir = fileparts(mfilename('fullpath'));
if nargin < 1, outfile = fullfile(rootdir, 'figures', 'loglik-sv.png'); end

%% data and posterior mean of the volatility path
read_data;                                           % YY, yall, tstart, tend, nlags
T = size(YY, 1);
XX = zeros(T, 2*nlags + 1);
for l = 1:nlags, XX(:, 2*l-1:2*l) = yall(tstart-l:tend-l, :); end
XX(:, end) = 1;
H = [];
for c = 1:4
    L = load(fullfile(rootdir, 'results', sprintf('chain_sv_%d.mat', c)), 'CH');
    H = [H L.CH.h]; %#ok<AGROW>
end
hbar = mean(H, 2);
w = exp(-hbar);
Pw = (XX'*(XX.*w))\(XX'*(YY.*w));                    % weighted least squares
Ew = YY - XX*Pw;
Om = (Ew.*w)'*Ew/T;                                  % Omega*

%% figure
navy = [20 43 141]/255; bg = [238 238 238]/255;      % background of the website (#eee)
set(groot, 'defaultAxesFontName', 'Helvetica', 'defaultTextFontName', 'Helvetica');
ga = linspace(-5, 5, 1401); gb = linspace(-5, 5, 1401);
[AG, BG] = meshgrid(ga, gb);
lev = -(100:-10:10);                                 % maximum - 100, ..., maximum - 10
f = figure('Units', 'inches', 'Position', [1 1 6.5 3.4], 'Color', bg, 'Visible', 'off');
set(f, 'DefaultAxesFontSize', get(groot, 'FactoryAxesFontSize'));   % layout independent of make_figures' defaults
d1 =BG.^2*Om(1,1) - 2*BG*Om(1,2) + Om(2,2);         % (A Omega A')_11, demand row (-beta, 1)
d2 = AG.^2*Om(1,1) - 2*AG*Om(1,2) + Om(2,2);         % (A Omega A')_22, supply row (-alpha, 1)
Z = T*log(abs(AG - BG)) - (T/2)*log(d1.*d2);
Zmax = -(T/2)*log(det(Om));                          % value on the ridge (constant along the curve)
Z(~isfinite(Z)) = NaN;
ax = axes(f); hold(ax, 'on');
contour(ax, AG, BG, Z - Zmax, 'LevelList', lev, 'LineColor', navy, 'LineWidth', 0.5);
plot(ax, [-5 5], [0 0], ':', 'Color', 'k', 'LineWidth', 0.5);
hL = Om(1,2)/Om(1,1);                                % vertical asymptote of the curve
for br = 1:2                                         % the two branches, not joined across the asymptote
    if br == 1, x = linspace(-5, hL - 1e-6, 20000); else, x = linspace(hL + 1e-6, 5, 20000); end
    y = (Om(2,2) - x*Om(1,2))./(Om(1,2) - x*Om(1,1));
    y(abs(y) > 5.05) = NaN;
    plot(ax, x, y, '-', 'Color', 'k', 'LineWidth', 1.6);
end
xlim(ax, [-5 5]); ylim(ax, [-5 5]);
set(ax, 'XTick', -5:5, 'YTick', [-5 0 5], 'Box', 'on', 'TickDir', 'in', 'TickLength', [0.01 0.01], ...
    'Color', bg, 'XColor', 'k', 'YColor', 'k', 'FontSize', 9, 'Layer', 'top', 'LineWidth', 0.5);
title(ax, 'Contours for log likelihood', 'FontWeight', 'normal', 'FontSize', 10, 'Color', 'k');
xlabel(ax, char(945), 'FontSize', 11, 'Interpreter', 'none');      % upright alpha and beta
ylabel(ax, char(946), 'FontSize', 11, 'Interpreter', 'none');
set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
print(f, outfile, '-dpng', '-r300');
close(f);
fprintf('Figure 1: Omega* = [%.4f %.4f; %.4f %.4f], residual correlation %.4f, asymptote alpha = %.4f\n', ...
    Om(1,1), Om(1,2), Om(2,1), Om(2,2), Om(1,2)/sqrt(Om(1,1)*Om(2,2)), hL);
end
