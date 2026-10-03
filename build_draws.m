function build_draws()
% BUILD_DRAWS  Pools the four chains of each model and computes, on 500 draws per
% chain, the unit and one-standard-deviation impulse responses, the forecast-error
% variance decomposition and the historical decomposition (contributions to growth).
% Saves results/draws.mat.
NUSE = 500; hmax = 21; Hf = 16;
rootdir = fileparts(mfilename('fullpath')); rdir = fullfile(rootdir,'results');
read_data;
T = size(YY,1); n = size(YY,2);
XX = []; for l = 1:nlags, XX = [XX yall(tstart-l:tend-l,:)]; end; XX = [XX ones(T,1)]; %#ok<AGROW>

R = struct();
for mdl = {'sv','nosv'}
    m = struct('a_all',[],'a',[],'B',[],'invd',[],'h',[],'phi',[],'sig2',[],'acc',[]);
    for c = 1:4
        L = load(fullfile(rdir, sprintf('chain_%s_%d.mat', mdl{1}, c)), 'CH'); CH = L.CH;
        m.a_all = [m.a_all CH.a]; m.phi = [m.phi CH.phi]; m.sig2 = [m.sig2 CH.sig2]; m.acc(c) = CH.acc;
        ix = round(linspace(1, size(CH.B,3), NUSE));     % B, D and h are stored every 10 draws
        m.a = [m.a CH.a(:,(ix-1)*10+1)]; m.B = cat(3,m.B,CH.B(:,:,ix));
        m.invd = [m.invd CH.invd(:,ix)]; m.h = [m.h CH.h(:,ix)];
    end
    nuse = size(m.a,2);
    m.irfU = zeros(n,n,hmax,nuse); m.irfS = zeros(n,n,hmax,nuse);
    m.fev  = zeros(n,n,Hf,nuse);   m.hdg  = zeros(n,n,T,nuse);
    tmp = zeros(n,n,hmax);
    for s = 1:nuse
        A = setA(m.a(:,s)); B = m.B(:,:,s); d = 1./m.invd(:,s);
        ma = structural_ma(A,B,nlags,n,T);
        m.irfU(:,:,:,s) = cumsum(ma(:,:,1:hmax),3);                 % unit shock, cumulated to levels
        for hh = 1:hmax, tmp(:,:,hh) = ma(:,:,hh)*diag(sqrt(d)); end
        m.irfS(:,:,:,s) = cumsum(tmp,3);                            % one-SD shock with e^{h_t} = 1
        contrib = zeros(n,n);
        for hh = 1:Hf, contrib = contrib + (ma(:,:,hh).^2).*(d'); m.fev(:,:,hh,s) = contrib./sum(contrib,2); end
        m.hdg(:,:,:,s) = compute_histdecomp(ma, YY*A' - XX*B');      % contributions to growth
    end
    R.(mdl{1}) = m;
end
R.YY = YY; R.time = time; R.T = T;
save(fullfile(rdir,'draws.mat'), 'R', '-v7.3');
fprintf('Saved results/draws.mat: %d draws of A per model, %d for IRF/FEVD/HD\n', size(R.sv.a_all,2), nuse);
end
