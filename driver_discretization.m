ns = [1000 2000 4000 8000];
gamma = 1e-3;
maxiter = 10000;
check = inf;
tol = 1e-10;

close all;
for n = [1000 2000 4000 8000]
    x = linspace(0,1,n)';
    mu = exp(-(x-0.2).^2*10^2) + 1*exp(-abs(x-0.4)*20)+0.01;
    mu = mu*(n-1)/sum(mu);
    nu = exp(-(x-0.6).^2*10^2)+0.01;
    nu = nu*(n-1)/sum(nu);
    cost = (x-x').^2;
    
    [pi,alpha,beta,errors,time,cgiters] = sinkhorn_newton(mu,nu,cost,gamma,maxiter,check,tol);

    figure(1)
    semilogy(cgiters,errors.mu); hold on; grid on;
    figure(2)
    semilogy(time,errors.mu); hold on; grid on; drawnow;
end

figure(1)
legend('N = 1000','N = 2000','N = 4000','N = 8000');
xlabel('errors over CG iterations');

figure(2)
legend('N = 1000','N = 2000','N = 4000','N = 8000');
xlabel('errors over run time');
