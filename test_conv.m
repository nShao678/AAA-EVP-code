close all
clear all
rng(1)

mMax = 8;
err = zeros(2,mMax);
ss = 100;
xx = exp(2*pi*1i*linspace(0,1,ss));
for m = 1:mMax
f = @(x) 1./((x-0.5).^m)./(x-2)./((x-3).^2);
[r,pols,~,zer] = aaa(f,xx);
pols = pols(abs(pols)<1);
err(1,m) = max(abs(pols-0.5));

err(2,m) = eps^(1/m);
end
latex(sym(err))



m = 8;
f = @(x) 1./((x-0.5).^m)./(x-2)./((x-3).^2);
rMax = 16;
base = 0.1;
err = zeros(2,rMax);
xx0 = exp(2*pi*1i*linspace(0,1,100));
for iterr = 1:rMax
    r = base^(iterr-1);
    xx = xx0*r+(1-r)/2;
[~,pols] = aaa(f,xx);
pols = pols(abs(pols)<1);
err(1,iterr) = max(abs(pols-0.5));
err(2,iterr) = eps^(1/m)*r;
end
latex(sym(err))


