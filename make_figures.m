% MAKE_FIGURES  The seven figures of the post (in figures/) and every number cited
% in the text and the table (printed). Reads results/draws.mat (build_draws.m).

rootdir = fileparts(mfilename('fullpath'));
load(fullfile(rootdir,'results','draws.mat'),'R');
out = fullfile(rootdir,'figures'); if ~isfolder(out), mkdir(out); end
time = R.time; T = R.T; sv = R.sv; ns = R.nosv; n = 2;
S = style();
save_png = @(f,name) exportgraphics(f, fullfile(out,name), 'Resolution',218, 'BackgroundColor',S.bg);
VN = {'Wages','Employment'}; SN = {'Demand','Supply'};

%% Figure 1: common volatility factor
V = exp(sv.h/2); qv = quantile(V,[0.05 0.16 0.5 0.84 0.95],2);
f = figure('Units','inches','Position',[1 1 7 3],'Color',S.bg,'Visible','off');
ax = axes(f); hold(ax,'on');
yl = [0 ceil(max(qv(:,5))*4)/4];
h90 = fill(ax,[time;flipud(time)],[qv(:,1);flipud(qv(:,5))],S.lblue2,'EdgeColor','none');
h68 = fill(ax,[time;flipud(time)],[qv(:,2);flipud(qv(:,4))],S.lblue,'EdgeColor','none');
hm  = plot(ax,time,qv(:,3),'-','Color',S.navy,'LineWidth',1.7);
yline(ax,1,'--','Color',S.zero,'LineWidth',0.8);
fmt_axes(ax,S); xlim(ax,[time(1) time(end)]); ylim(ax,yl); xticks(ax,1970:10:2020);
ylabel(ax,'$e^{h_t/2}$','Interpreter','latex','FontSize',10);
lg = legend(ax,[hm h68 h90],{'Median','68% band','90% band'},'Location','northoutside', ...
    'Orientation','horizontal','Box','off','FontSize',8.5); %#ok<NASGU>
save_png(f,'volatility-factor.png');

%% Figure 2: posterior histograms of the elasticities and the prior
nm = {'Demand slope \beta','Supply slope \alpha'};
f = figure('Units','inches','Position',[1 1 7 2.9],'Color',S.bg,'Visible','off');
t = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
cA = [-0.6 0.6]; sA = 0.6; nu = 3;
grids = {linspace(-4,0,600), linspace(0,1.2,600)};
for ia = 1:2
    ax = nexttile(t); hold(ax,'on'); x = grids{ia};
    sg = sign(cA(ia));
    pr = tpdf((x-cA(ia))/sA,nu)/sA; pr = pr / (1 - tcdf(-sg*cA(ia)/sA,nu));   % sign-truncated Student-t
    e = linspace(x(1), x(end), 81);                          % same bins for both models
    hb = histogram(ax, ns.a_all(ia,:), 'BinEdges',e, 'Normalization','pdf', ...
        'FaceColor',S.gband,'FaceAlpha',0.9,'EdgeColor','none');
    hs = histogram(ax, sv.a_all(ia,:), 'BinEdges',e, 'Normalization','pdf', ...
        'FaceColor',S.lblue,'FaceAlpha',0.6,'EdgeColor','none');
    hp = plot(ax,x,pr,'-','Color',[0.2 0.2 0.2],'LineWidth',1.0);
    kb = hb.Values; ks = hs.Values;
    fmt_axes(ax,S); xlim(ax,[x(1) x(end)]); ylim(ax,[0 max([kb ks])*1.08]);
    title(ax,nm{ia},'FontWeight','normal','FontSize',9.5,'Color',S.ax);
    if ia==1, ylabel(ax,'Density'); end
end
lg = legend([hs hb hp],{'Posterior, SV','Posterior, no SV','Prior'}, ...
    'Orientation','horizontal','Box','off','FontSize',8.5); lg.Layout.Tile = 'south';
save_png(f,'elasticity-posteriors.png');

%% Figures 3 and 4: historical decomposition, contributions to 8-quarter growth (SV model)
W = 8;
win = @(m,i) movsum(squeeze(m.hdg(i,:,:,:)), [W-1 0], 2, 'Endpoints','fill');   % n x T x draws
names = {'historical-decomposition-wages.png','historical-decomposition.png'};
for i = 1:2
    f = single_hd_figure(time, win(sv,i), {'Demand shocks','Supply shocks'}, 'Percentage points', S);
    save_png(f, names{i});
end

%% Figure 5: responses to a unit structural shock
hmax = size(sv.irfU,3); HO = (0:hmax-1)';
f = figure('Units','inches','Position',[1 1 7 5],'Color',S.bg,'Visible','off');
t = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for i=1:n, for j=1:n
    ax = nexttile(t); hold(ax,'on');
    qs = quantile(squeeze(sv.irfU(i,j,:,:)),[0.16 0.5 0.84],2);
    qb = quantile(squeeze(ns.irfU(i,j,:,:)),[0.16 0.5 0.84],2);
    hh = draw_compare(ax, HO, qs, qb, S);
    fmt_axes(ax,S); xlim(ax,[0 hmax-1]); xticks(ax,0:4:20);
    title(ax,sprintf('%s, %s shock',VN{i},lower(SN{j})),'FontWeight','normal','FontSize',9.5,'Color',S.ax);
    if i==n, xlabel(ax,'Quarters'); end
    if j==1, ylabel(ax,'Percent'); end
end, end
lg = legend(hh,{'SV, median','SV, 68% band','No SV, median','No SV, 68% band'}, ...
    'Orientation','horizontal','Box','off','FontSize',8.5); lg.Layout.Tile = 'south';
save_png(f,'unit-irfs.png');

%% Figure 6: one-standard-deviation responses of employment, 2009 and 1995
[~,ihi] = min(abs(time-2009)); [~,ilo] = min(abs(time-1995));
ehi = exp(sv.h(ihi,:)/2); elo = exp(sv.h(ilo,:)/2);
f = figure('Units','inches','Position',[1 1 7 2.9],'Color',S.bg,'Visible','off');
t = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for j = 1:2
    ax = nexttile(t); hold(ax,'on');
    shp = squeeze(sv.irfS(2,j,:,:));
    qh = quantile(shp.*ehi,[0.16 0.5 0.84],2); ql = quantile(shp.*elo,[0.16 0.5 0.84],2);
    a1 = fill(ax,[HO;flipud(HO)],[qh(:,1);flipud(qh(:,3))],S.lblue,'EdgeColor','none','FaceAlpha',0.6);
    a2 = fill(ax,[HO;flipud(HO)],[ql(:,1);flipud(ql(:,3))],S.gband,'EdgeColor','none','FaceAlpha',0.8);
    l1 = plot(ax,HO,qh(:,2),'-','Color',S.navy,'LineWidth',1.5);
    l2 = plot(ax,HO,ql(:,2),'-','Color',[0.2 0.2 0.2],'LineWidth',1.5);
    yline(ax,0,'--','Color',S.zero,'LineWidth',0.8);
    fmt_axes(ax,S); xlim(ax,[0 hmax-1]); xticks(ax,0:4:20); xlabel(ax,'Quarters');
    title(ax,sprintf('Employment, %s shock',lower(SN{j})),'FontWeight','normal','FontSize',9.5,'Color',S.ax);
    if j==1, ylabel(ax,'Percent'); end
end
lg = legend([l1 a1 l2 a2],{'2009, median','2009, 68% band','1995, median','1995, 68% band'}, ...
    'Orientation','horizontal','Box','off','FontSize',8.5); lg.Layout.Tile = 'south';
save_png(f,'one-sd-employment.png');

%% Figure 7: ACF and PACF of the squared structural shocks, draw by draw
read_data;
XX = []; for l = 1:nlags, XX = [XX yall(tstart-l:tend-l,:)]; end; XX = [XX ones(T,1)]; %#ok<AGROW>
K = 20; cb = 1.96/sqrt(T); nd = size(sv.a,2);
AC = zeros(K,2,2,nd); PA = zeros(K,2,2,nd); LBp = zeros(2,2,nd);      % (lag, shock, model, draw)
mods = {sv, ns};
for mm = 1:2
    m = mods{mm};
    for s = 1:nd
        A = setA(m.a(:,s)); Bm = m.B(:,:,s); d = 1./m.invd(:,s);
        if mm == 1, sc = sqrt(d'.*exp(m.h(:,s))); else, sc = sqrt(d'); end
        E2 = ((YY*A' - XX*Bm')./sc).^2;
        for i = 1:2
            AC(:,i,mm,s) = myacf(E2(:,i),K); PA(:,i,mm,s) = mypacf(E2(:,i),K);
            [~,LBp(i,mm,s)] = ljungbox(E2(:,i),K);
        end
    end
end
qAC = quantile(AC,[0.16 0.5 0.84],4); qPA = quantile(PA,[0.16 0.5 0.84],4);
f = figure('Units','inches','Position',[1 1 7 4.8],'Color',S.bg,'Visible','off');
t = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
lag = (1:K)'; off = 0.18; gcol = [0.45 0.45 0.45]; hmk = gobjects(1,2);
for i = 1:2
    for c = 1:2
        ax = nexttile(t); hold(ax,'on');
        if c==1, Q = qAC; ttl = 'ACF'; else, Q = qPA; ttl = 'PACF'; end
        % +-1.96/sqrt(T) band
        patch(ax,[0.3 K+0.7 K+0.7 0.3],[-cb -cb cb cb],[0.86 0.86 0.86],'EdgeColor','none','HandleVisibility','off');
        yline(ax,0,'-','Color',[0.55 0.55 0.55],'LineWidth',0.6);
        for mm = [2 1]                                           % no SV on the left, SV on the right
            xj = lag + (2*(mm==1)-1)*off;
            lo = squeeze(Q(:,i,mm,1)); md = squeeze(Q(:,i,mm,2)); hi = squeeze(Q(:,i,mm,3));
            if mm == 1, col = S.navy; fc = S.navy; else, col = gcol; fc = S.bg; end
            plot(ax,[xj xj]',[lo hi]','-','Color',col,'LineWidth',1.0);      % 68% interval
            hmk(mm) = plot(ax,xj,md,'o','MarkerSize',4,'Color',col,'MarkerFaceColor',fc,'LineWidth',1.0);
        end
        fmt_axes(ax,S); xlim(ax,[0.3 K+0.7]); ylim(ax,[-0.3 0.45]); xticks(ax,[1 4:4:20]);
        title(ax,sprintf('%s shock squared: %s', SN{i}, ttl),'FontWeight','normal','FontSize',9.5,'Color',S.ax);
        if i==2, xlabel(ax,'Lag'); end
    end
end
lg = legend(hmk([2 1]),{'No SV, median and 68% interval','SV, median and 68% interval'}, ...
    'Orientation','horizontal','Box','off','FontSize',8.5);
lg.Layout.Tile = 'south';
save_png(f,'acf-squared-shocks.png');
r1 = squeeze(AC(1,:,:,:));                                     % shock x model x draw

%% Numbers cited in the text
fprintf('\n===== Sample %.2f-%.2f, T=%d =====\n', time(1), time(end), T);
fprintf('phi mean %.3f [%.3f %.3f], sigma2 mean %.4f\n', mean(sv.phi), quantile(sv.phi,.05), quantile(sv.phi,.95), mean(sv.sig2));
md = qv(:,3);
for r = {[1970 1978],[1979 1984],[1985 2007],[2007 2012],[2012 2020]}
    w = time>=r{1}(1) & time<r{1}(2);
    [mx,imx] = max(md.*w./w); [mn,imn] = min(md + 1e9*(~w)); tt = time(w);
    fprintf('  %d-%d: max %.2f at %.2f, min %.2f at %.2f\n', r{1}, mx, time(imx), mn, time(imn));
end
fprintf('  first quarter below 1 after 1984: %.2f\n', time(find(time>1984 & md<1,1)));
for f2 = {'nosv','sv'}
    m = R.(f2{1});
    for ia = 1:2
        x = m.a_all(ia,:);
        fprintf('%-5s %s median %6.2f sd %5.2f 68%% [%6.2f %6.2f] 95%% [%6.2f %6.2f]\n', f2{1}, nm{ia}, ...
            median(x), std(x), quantile(x,[.16 .84]), quantile(x,[.025 .975]));
    end
end
Hf = size(sv.fev,3);
for i=1:2
    for k=1:2
        s = squeeze(sv.fev(i,k,Hf,:)); b = squeeze(ns.fev(i,k,Hf,:));
        fprintf('FEVD h=%d %s share in %s: SV %.2f [%.2f %.2f]  noSV %.2f [%.2f %.2f]\n', Hf, SN{k}, VN{i}, ...
            median(s), quantile(s,[.16 .84]), median(b), quantile(b,[.16 .84]));
    end
end
bwF = @(m) mean(quantile(m.fev(:,2,Hf,:),0.84,4)-quantile(m.fev(:,2,Hf,:),0.16,4),'all');
fprintf('FEVD supply-share width ratio SV/noSV: %.2f\n', bwF(sv)/bwF(ns));
bw = @(X) squeeze(quantile(X,0.84,3)-quantile(X,0.16,3));
for i=1:2
    Ws = bw(win(sv,i)); Wb = bw(win(ns,i));
    fprintf('HD (8 quarters) %s width ratio SV/noSV %.2f (demand %.2f, supply %.2f)\n', VN{i}, ...
        mean(Ws(:),'omitnan')/mean(Wb(:),'omitnan'), mean(Ws(1,:),'omitnan')/mean(Wb(1,:),'omitnan'), ...
        mean(Ws(2,:),'omitnan')/mean(Wb(2,:),'omitnan'));
end
% troughs and peaks of the median contributions (SV model) around each episode,
% and the largest gap between the SV and homoskedastic medians
episodes = [1974 1977; 1980 1985; 1990 1994; 2001 2005; 2008 2012; 1998 2001];
for i = 1:2
    Ms = median(win(sv,i),3); Mb = median(win(ns,i),3);
    fprintf('%s: max |SV - no SV| median gap: demand %.2f, supply %.2f\n', VN{i}, ...
        max(abs(Ms(1,:)-Mb(1,:)),[],'omitnan'), max(abs(Ms(2,:)-Mb(2,:)),[],'omitnan'));
    for e = 1:size(episodes,1)
        w = time >= episodes(e,1) & time < episodes(e,2);
        for kk = 1:2
            v = Ms(kk,:); v(~w) = NaN;
            [mn,jmn] = min(v); [mx,jmx] = max(v);
            fprintf('   %d-%d %-6s min %6.2f (%.2f)  max %6.2f (%.2f)\n', episodes(e,:), SN{kk}, mn, time(jmn), mx, time(jmx));
        end
    end
end
for i=1:2, for j=1:2
    fprintf('unit IRF %s<-%s impact SV %.2f noSV %.2f | h=20 SV %.2f noSV %.2f\n', VN{i}, SN{j}, ...
        median(sv.irfU(i,j,1,:)), median(ns.irfU(i,j,1,:)), median(sv.irfU(i,j,end,:)), median(ns.irfU(i,j,end,:)));
end, end
fprintf('one-SD: exp(h/2) 2009 %.2f 1995 %.2f ratio %.2f\n', median(ehi), median(elo), median(ehi)/median(elo));
mn = {'SV','noSV'};
for i=1:2, for mm=1:2
    x = squeeze(r1(i,mm,:));
    fprintf('%s squared %-4s: r1 median %.2f [16-84: %.2f %.2f] P(r1>cb)=%.2f | LB(%d) p median %.3f\n', SN{i}, mn{mm}, ...
        median(x), quantile(x,[.16 .84]), mean(x>cb), K, median(squeeze(LBp(i,mm,:))));
end, end
fprintf('Figures saved in %s\n', out);

%% Local functions
function f = single_hd_figure(time, Hs, ttl, ylab, S)
% SV model only: median with 68% and 90% bands
f = figure('Units','inches','Position',[1 1 7 3.1],'Color',S.bg,'Visible','off');
t = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
Q = squeeze(quantile(Hs,[0.05 0.16 0.5 0.84 0.95],3));       % shock x T x 5
ok = ~isnan(Q(1,:,3))'; x = time(ok);
lims = [min(Q(:)) max(Q(:))]; lims = lims + [-1 1]*0.05*diff(lims);
for k = 1:2
    ax = nexttile(t); hold(ax,'on'); q = squeeze(Q(k,ok,:));
    h90 = fill(ax,[x;flipud(x)],[q(:,1);flipud(q(:,5))],S.lblue2,'EdgeColor','none');
    h68 = fill(ax,[x;flipud(x)],[q(:,2);flipud(q(:,4))],S.lblue,'EdgeColor','none');
    hm  = plot(ax,x,q(:,3),'-','Color',S.navy,'LineWidth',1.6);
    yline(ax,0,'--','Color',S.zero,'LineWidth',0.8);
    fmt_axes(ax,S); xlim(ax,[time(1) time(end)]); ylim(ax,lims); xticks(ax,1970:10:2020);
    title(ax, ttl{k}, 'FontWeight','normal','FontSize',9.5,'Color',S.ax);
    if k == 1, ylabel(ax, ylab); end
end
lg = legend([hm h68 h90], {'Median','68% band','90% band'}, ...
    'Orientation','horizontal','Box','off','FontSize',8.5); lg.Layout.Tile = 'south';
end

function S = style()
S.navy = [0.03 0.17 0.45]; S.lblue = [0.58 0.72 0.88]; S.lblue2 = [0.82 0.88 0.95];
S.grey = [0.45 0.45 0.45]; S.gband = [0.80 0.80 0.80]; S.ax = [0.15 0.15 0.15]; S.zero = [0.25 0.25 0.25];
S.bg = [238 238 238]/255;                        % background of the website (#eee)
set(groot,'defaultAxesFontName','Helvetica','defaultTextFontName','Helvetica', ...
    'defaultLegendFontName','Helvetica','defaultAxesFontSize',8.5);
end

function fmt_axes(ax,S)
grid(ax,'on');
set(ax,'Box','off','Color',S.bg,'FontSize',8.5,'GridAlpha',0.15,'XColor',S.ax,'YColor',S.ax,'Layer','top');
end

function hh = draw_compare(ax, x, qs, qb, S)
hb1 = fill(ax,[x;flipud(x)],[qb(:,1);flipud(qb(:,3))],S.gband,'EdgeColor','none');
hs  = fill(ax,[x;flipud(x)],[qs(:,1);flipud(qs(:,3))],S.lblue,'EdgeColor','none','FaceAlpha',0.6);
hbm = plot(ax,x,qb(:,2),'--','Color','k','LineWidth',1.2);
hsm = plot(ax,x,qs(:,2),'-','Color',S.navy,'LineWidth',1.6);
yline(ax,0,'--','Color',S.zero,'LineWidth',0.8);
hh = [hsm hs hbm hb1];
end
