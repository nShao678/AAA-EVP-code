rng(1)
n = 100;
r = 2;
A1 = orth(randn(n,r));
A2 = orth(randn(n,r))';

sampleMax = 10000000;
bMax = 2;
err1 = zeros(sampleMax,bMax);
err2 = zeros(sampleMax,bMax);
err3 = zeros(sampleMax,bMax);
D1 = diag([1,1]);
D2 = diag([1,1e-2]);
D3 = diag([1,1e-4]);
for sample = 1:sampleMax
%     if mod(sample,10000)==0
%         sample
%     end
    OmegaL = randn(n,bMax);
    OmegaR = randn(n,bMax);
    S1 = (OmegaL'*A1)*D1*(A2*OmegaR);
    S2 = (OmegaL'*A1)*D2*(A2*OmegaR);
    S3 = (OmegaL'*A1)*D3*(A2*OmegaR);
    for b = 1:bMax
        err1(sample,b) = norm(S1(1:b,1:b),'fro');
        err2(sample,b) = norm(S2(1:b,1:b),'fro');
        err3(sample,b) = norm(S3(1:b,1:b),'fro');
    end
end
%%
err = [err1,err2,err3];
iiMax = 10;
base = 0.6;
data = zeros(bMax*3,iiMax);
for ii = 1:iiMax
    eps = base^(ii);
    data(:,ii) = sum(err<eps)/sampleMax;
end
figure
hold on

plot(base.^(1:iiMax),data(1,:),'r-o','LineWidth',2,'DisplayName','$b=1,\sigma_{2}=1$')
plot(base.^(1:iiMax),data(3,:),'r-x','LineWidth',2,'DisplayName','$b=1,\sigma_{2}=10^{-2}$')
plot(base.^(1:iiMax),data(5,:),'r-+','LineWidth',2,'DisplayName','$b=1,\sigma_{2}=10^{-4}$')


plot(base.^(1:iiMax),data(2,:),'b-o','LineWidth',2,'DisplayName','$b=2,\sigma_{2}=1$')
plot(base.^(1:iiMax),data(4,:),'b-x','LineWidth',2,'DisplayName','$b=2,\sigma_{2}=10^{-2}$')
plot(base.^(1:iiMax),data(6,:),'b-+','LineWidth',2,'DisplayName','$b=2,\sigma_{2}=10^{-4}$')

plot(base.^(1:iiMax),base.^((1:iiMax))/2,'k--d','LineWidth',2,'DisplayName','$y\propto\epsilon^{-1}$')
plot(base.^(1:iiMax),base.^(2*(1:iiMax))*2,'k--d','LineWidth',2,'DisplayName','$y\propto\epsilon^{-2}$')
plot(base.^(1:iiMax),base.^(3*(1:iiMax)),'k--d','LineWidth',2,'DisplayName','$y\propto\epsilon^{-3}$')


hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Interpreter','latex','FontSize',12,'location','southeast','box','off')
set(gcf, 'Color', 'w');
export_fig(['fig/gaussianr.pdf'])
export_fig(['fig/gaussianr.eps'])


%% 

rng(1)
n = 100;
r = 2;
A1 = orth(randn(n,r)+1i*randn(n,r));
A2 = orth(randn(n,r)+1i*randn(n,r))';

sampleMax = 10000000;
bMax = 2;
err1 = zeros(sampleMax,bMax);
err2 = zeros(sampleMax,bMax);
err3 = zeros(sampleMax,bMax);
D1 = diag([1,1]);
D2 = diag([1,1e-2]);
D3 = diag([1,1e-4]);
for sample = 1:sampleMax
    OmegaL = randn(n,bMax)+1i*randn(n,bMax);
    OmegaR = randn(n,bMax)+1i*randn(n,bMax);
    S1 = (OmegaL'*A1)*D1*(A2*OmegaR);
    S2 = (OmegaL'*A1)*D2*(A2*OmegaR);
    S3 = (OmegaL'*A1)*D3*(A2*OmegaR);
    for b = 1:bMax
        err1(sample,b) = norm(S1(1:b,1:b),'fro');
        err2(sample,b) = norm(S2(1:b,1:b),'fro');
        err3(sample,b) = norm(S3(1:b,1:b),'fro');
    end
end
%%
err = [err1,err2,err3];
iiMax = 10;
base = 0.85;
data = zeros(bMax*3,iiMax);
for ii = 1:iiMax
    eps = base^(ii);
    data(:,ii) = sum(err<eps)/sampleMax;
end
figure
hold on

plot(base.^(1:iiMax),data(1,:),'r-o','LineWidth',2,'DisplayName','$b=1,\sigma_{2}=1$')
plot(base.^(1:iiMax),data(3,:),'r-x','LineWidth',2,'DisplayName','$b=1,\sigma_{2}=10^{-2}$')
plot(base.^(1:iiMax),data(5,:),'r-+','LineWidth',2,'DisplayName','$b=1,\sigma_{2}=10^{-4}$')


plot(base.^(1:iiMax),data(2,:),'b-o','LineWidth',2,'DisplayName','$b=2,\sigma_{2}=1$')
plot(base.^(1:iiMax),data(4,:),'b-x','LineWidth',2,'DisplayName','$b=2,\sigma_{2}=10^{-2}$')
plot(base.^(1:iiMax),data(6,:),'b-+','LineWidth',2,'DisplayName','$b=2,\sigma_{2}=10^{-4}$')

plot(base.^(1:iiMax),base.^(2*(1:iiMax))/10,'k--d','LineWidth',2,'DisplayName','$y\propto\epsilon^{-2}$')
plot(base.^(1:iiMax),base.^(4*(1:iiMax))/50,'k--d','LineWidth',2,'DisplayName','$y\propto\epsilon^{-4}$')
plot(base.^(1:iiMax),base.^(6*(1:iiMax))/1000,'k--d','LineWidth',2,'DisplayName','$y\propto\epsilon^{-6}$')


hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Interpreter','latex','FontSize',12,'location','southeast','box','off')
set(gcf, 'Color', 'w');
export_fig(['fig/gaussianc.pdf'])
export_fig(['fig/gaussianc.eps'])