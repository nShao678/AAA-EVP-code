close all
rng(1)
sample = 50;
center = 0.7+0.6*1i;
rad = 0.5;
NMax = 5;
bMax = 3;
para.fval = 1;
err = zeros(NMax,bMax);
[coef,~,T] = butterfly;
n = size(coef{1},1);
OmegaL0 = randn(n,bMax);
OmegaR0 = randn(n,bMax);
zzSet = cell(NMax,bMax);
fzSet = cell(NMax,bMax);

% compute reference sol
region = center+rad*([1+1i,1-1i,-1-1i,-1+1i]);
eval_p0 = polyeig(coef{:});
eval_p = eval_p0(inpoly(eval_p0,region));

%% plot region
figure
hold on

plot(eval_p,'ro','MarkerSize',8,'LineWidth',1)
plot(eval_p0,'bx','MarkerSize',8,'LineWidth',1)
plot(contour(center,rad,sample),'k.')

hold off
set(gcf, 'Color', 'w');
export_fig(['fig/butterfly1.pdf'])
export_fig(['fig/butterfly1.eps'])

figure
hold on

plot(eval_p,'ro','MarkerSize',8,'LineWidth',1)
plot(eval_p0,'bx','MarkerSize',8,'LineWidth',1)

rads = rad/3;
for ii = 1:3
    for jj = 1:3
        centers = center-rad*(1+1i)+2*ii*rads-rads+1i*(2*jj*rads-rads);        
        plot(contour(centers,rads,sample),'k.')
    end
end
hold off
set(gcf, 'Color', 'w');
axis([0.1,1.3,0,1.2])
export_fig(['fig/butterfly2.pdf'])
export_fig(['fig/butterfly2.eps'])


%% aaa_evp

for b = 1:bMax
    for N = 1:NMax
        OmegaL = OmegaL0(:,1:b);
        OmegaR = OmegaR0(:,1:b);
        fun = @(z) OmegaL'*(T(z)\OmegaR);
        [pol,zz,ff] = aaa_evp(fun,center,rad,sample,b,N,para);
        if para.fval == 1
            fz = zeros(size(ff));
            for ii = 1:size(ff)
                fz(ii) = norm(ff{ii}-fun(zz(ii)));
            end
            zzSet{N,b} = zz;
            fzSet{N,b} = fz;

        end
        err(N,b) = cmpvec(eval_p,pol);
    end
end
%% plot
figure
t = tiledlayout(bMax, NMax, 'TileSpacing', 'compact', 'Padding', 'loose'); % balance

for b = 1:bMax
    for N = 1:NMax
        nexttile;
        hold on
        plot_fx(zzSet{N,b}, fzSet{N,b}, center, rad);
        plot(eval_p, 'wx');
        plot(eval_p, 'ko');
        hold off
        % title(['N=',num2str(N),' and b=',num2str(b)],'FontSize',18)
        title(['dist = ', num2str(err(N,b), '%0.0e')], 'FontSize', 16);
        set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    end
end

colormap jet;
cb = colorbar;
cb.Layout.Tile = 'south';
cb.FontSize = 16;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.85]); % large & readable
export_fig(['fig/butterfly.pdf'])
export_fig(['fig/butterfly.eps'])


%% Fix degree
clear all
rng(1)
sample = 50;
center = 0.5+0.5*1i;
rad = 0.2;
bMax = 5;
para.fval = 1;
err = zeros(2,bMax);
deg = zeros(2,bMax);
[coef,~,T] = butterfly;
n = size(coef{1},1);
OmegaL0 = randn(n,bMax);
OmegaR0 = randn(n,bMax);
zzSet = cell(2,bMax);
fzSet = cell(2,bMax);
% compute reference sol
region = center+rad*([1+1i,1-1i,-1-1i,-1+1i]);
eval_p0 = polyeig(coef{:});
eval_p = eval_p0(inpoly(eval_p0,region));
% figure
% hold on
% 
% plot(eval_p,'ro','MarkerSize',8,'LineWidth',1)
% plot(eval_p0,'bx','MarkerSize',8,'LineWidth',1)
% plot(contour(center,rad,sample),'k.')
% 
% hold off
% set(gcf, 'Color', 'w');

for iter = 1:2
    switch iter
        case 1
        para.deg = 70;
        case 2
        para.deg = 51;
    end
    for b = 1:bMax
        OmegaL = OmegaL0(:,1:b);
        OmegaR = OmegaR0(:,1:b);
        fun = @(z) OmegaL'*(T(z)\OmegaR);
        [pol,zz,ff,deg(iter,b)] = aaa_evp(fun,center,rad,sample,b,1,para);
        if para.fval == 1
            fz = zeros(size(ff));
            for ii = 1:size(ff)
                fz(ii) = norm(ff{ii}-fun(zz(ii)));
            end
            zzSet{iter,b} = zz;
            fzSet{iter,b} = fz;
        end
        err(iter,b) = cmpvec(eval_p,pol);
    end
end
%% plot
figure
t = tiledlayout(2,bMax, 'TileSpacing', 'compact', 'Padding', 'loose'); % balance


for iter = 1:2
    for b = 1:bMax
        nexttile;
        hold on
        plot_fx(zzSet{iter,b}, fzSet{iter,b}, center, rad);
        plot(eval_p, 'wx');
        plot(eval_p, 'ko');
        hold off
        % title(['N=',num2str(N),' and b=',num2str(b)],'FontSize',18)
        title(['dist = ', num2str(err(iter,b), '%0.0e'), ' degree=',num2str(deg(iter,b))], 'FontSize', 16);
        set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    end
end

colormap jet;
cb = colorbar;
cb.Layout.Tile = 'south';
cb.FontSize = 16;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.85]); % large & readable

export_fig(['fig/butterflyDeg.pdf'])
export_fig(['fig/butterflyDeg.eps'])

function [pol,zz,ff,deg] = aaa_evp(fun,center,rad,sample,b,N,para)
% deg only works when N=1
pol = [];
zz = [];
ff = [];
rads = rad/N;
if para.fval==1
    numplot = 100;
    zz0 = linspace(-rads,rads,numplot);
    zz0 = reshape(zz0'*ones(1,numplot)*1i+ones(numplot,1)*zz0,[],1);
end
for ii = 1:N
    for jj = 1:N
        centers = center-rad*(1+1i)+2*ii*rads-rads+1i*(2*jj*rads-rads);
        Z = contour(centers,rads,sample);
        samples = length(Z);
        F = zeros(samples,b^2);
        for kk = 1:samples
            F(kk,:) = reshape(fun(Z(kk)),1,b^2);
        end
        nF = norm(F,'fro');
        F = F/nF;
        if isfield(para,'deg')
            [r,pols] = sv_aaa_mod(F,Z,'mmax',para.deg);
        else
            [r,pols] = sv_aaa_mod(F,Z);
        end
        deg = length(pols);
        pol = [pol;pols(inpoly(pols,Z))];
        if para.fval==1
            zzs = zz0+centers;
            ffs = cell(length(zzs),1);
            for kk = 1:length(zzs)
                ffs{kk} = nF*reshape(r(zzs(kk)),b,b);
            end
            zz = [zz;zzs];
            ff = [ff;ffs];
        end
    end
end


end

function idx = inpoly(z,w)
idx = inpolygon(real(z),imag(z),real(w),imag(w));
end
function zz = contour(center,rad,sample)

zz0 = linspace(-rad,rad-2*rad/sample,sample);
zz = [zz0+center+rad*1i,-1i*zz0+center+rad,-zz0+center-rad*1i,1i*zz0+center-rad];
end
function err = cmpvec(x,y)
err = 0;
if isempty(y)
    err = inf;
else
    for ii = 1:length(x)
        err = max(err,min(abs(x(ii)-y)));
    end
end
end

function plot_fx(x, fx,center,rad)

N = 300; % resolution
xx = real(center)+linspace(-rad, rad, N);
yy = imag(center)+linspace(-rad, rad, N);
[X, Y] = meshgrid(xx,yy);
Zc = X + 1i*Y;
% mask = abs(Zc) <= 1; % inside the unit disk

% Interpolate log|f(x)| onto the grid
vals = log10(abs(fx(:)) + eps); % avoid log(0)
F = scatteredInterpolant(real(x(:)), imag(x(:)), vals, 'natural', 'none');
FZ = F(real(Zc), imag(Zc));
% FZ(~mask) = NaN; % outside disk set to NaN

% Use imagesc instead of pcolor (no missing row/col)
imagesc(xx,yy, FZ);
set(gca, 'YDir', 'normal'); % keep y increasing upward
clim([-14,0])
axis equal tight off;
end