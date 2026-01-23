rng(1)
n = 100;
e = ones(n,1);
table = zeros(2,7);
for iter = 1:7
scale = 10^iter*100;
z = (1:n)';
z(1) = scale;
z(2) = scale+1;
z = z/scale;
w = [1,0.1,1e-12*ones(1,n-2)];
w = w/norm(w);
A = [0,w;e,diag(z)];
B = eye(n+1);
B(1,1) = 0;
[V,D] = eig(A,B,'vector');
idx = 3;
zz = 1./(D(idx)-z);
table(1,iter) = max(abs(D(idx)-z(1:2)));
table(2,iter) = abs(w*zz)/norm(zz);
end

latex(sym(table))
mMax = 3;
err = zeros(2,mMax+1);
for m = 0:mMax
omega = exp(2^(-m)*pi*1i);
z = omega.^(1:2^(m+1));
w = z.^(2^m+1);
e = ones(2^(m+1),1);
A = [0,w;e,diag(z)];
B = eye(2^(m+1)+1);
B(1) = 0;
eval = eig(A,B);
err(1,m+1) = min(abs(eval(~isinf(eval))));
err(2,m+1) = eps^(2^(-m));

end
latex(sym(err))