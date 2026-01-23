close all
rng(1)
eval = 2;
A = eval*eye(10);
A(1,2) = 1;
A(2,3) = 1;
A(4,5) = 1;
A(5,6) = 1;
A(7,8) = 1;
n = size(A,1);
b = 6;
iterMax = 15;
% OmegaL = randn(n,b);
% OmegaR = randn(n,b);

OmegaL = randn(n,b)+1i*randn(n,b);
OmegaR = randn(n,b)+1i*randn(n,b);
sval = zeros(b,iterMax);

for iter = 1:iterMax
    rho = eval+0.5^(iter);
    F = OmegaL'*((A-rho*eye(n))\OmegaR);
    sval(:,iter) = svd(F);
end
xx = 0.5.^(1:iterMax);
figure
hold on
for iter = 1:b
    plot(xx,sval(iter,:),'-o','LineWidth',3,'DisplayName',['$\sigma_{',num2str(iter),'}$'])
end
xx = 0.5.^(1:iterMax);
plot(xx,xx.^(-3),'k--x','LineWidth',2,'DisplayName','$y=\epsilon^{-3}$')
plot(xx,xx.^(-2),'b--x','LineWidth',2,'DisplayName','$y=\epsilon^{-2}$')
plot(xx,xx.^(-1),'r--x','LineWidth',2,'DisplayName','$y=\epsilon^{-1}$')
hold off
legend('FontSize',14,'Interpreter','latex','Location','northeast','Box','off','NumColumns',3);
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gcf, 'Color', 'w');
export_fig(['fig/multi.pdf'])
export_fig(['fig/multi.eps'])
