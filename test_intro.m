rng(1)
n = 9;
sample = 50;
A = diag(linspace(0.1,0.9,n));
z = exp(2*pi*1i*(1:sample)/sample)'/2+0.5;
omegaL = randn(n,1)+1i*randn(n,1);
omegaR = randn(n,1)+1i*randn(n,1);
f = zeros(sample,1);
for ii = 1:sample
    f(ii) = omegaL'*((A-z(ii)*eye(n))\omegaR);
end
[~,pol] = aaa(f,z);
sort(real(pol))


sample = 50;
omegaL = randn(n,2)+1i*randn(n,2);
omegaR = randn(n,2)+1i*randn(n,2);
f = zeros(sample,4);
pol = [];
for iter = 1:3
    center = (2*iter-1)/6;
z = exp(2*pi*1i*(1:sample)/sample)'/6+center;
for ii = 1:sample
    f(ii,:) = reshape(omegaL'*((A-z(ii)*eye(n))\omegaR),[],4);
end
[~,polt] = sv_aaa_mod(f,z);
pol = [pol;polt(abs(polt-center)<=1/6)];
end
sort(real(pol))