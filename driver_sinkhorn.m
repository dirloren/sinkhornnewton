n = 20;

% produce two input images of unit mass
[X,Y] = meshgrid(linspace(0,1,n));
img1 = exp(-((X-1/3).^2+(Y-1/3).^2)*6^2)+0.1;
img1 = img1/sum(img1(:));
img2 = exp(-((X-2/3).^2+(Y-2/3).^2)*3^2)+0.1;
img2 = img2/sum(img2(:));

% Wasserstein-p cost
p = 2;
cost = abs( repmat(X(:),1,length(X(:))) - repmat(X(:)',length(X(:)),1) ).^p + abs( repmat(Y(:),1,length(Y(:))) - repmat(Y(:)',length(Y(:)),1) ).^p;

% algorithm parameters
gamma = 0.001;
maxiter = 10000;
check = inf;
tol = 1e-13;

% solve once to get optimal plan and cost
[piopt] = sinkhorn_newton_primal(img1(:),img2(:),cost,gamma,maxiter,check,tol);
optcost = sum(sum(cost.*piopt));

% solve optimal transport problem
[pi,errors2,time2,cgiters] = sinkhorn_newton_primal(img1(:),img2(:),cost,gamma,maxiter,check,tol,'optplan',piopt,'optcost',optcost);
[~,v,w,errors,time] = sinkhorn(img1(:),img2(:),cost,gamma,maxiter,check,tol,'optplan',piopt,'optcost',optcost);

% plot
close all
figure(1); 
iters = 0:length(errors.mu)-1;
semilogy(iters,errors.mu,  '--',...
         iters,errors.cost,'--',...
         iters,errors.pi,  '--'); hold on;
ax = gca; ax.ColorOrderIndex = 1;
semilogy(cgiters,errors2.mu,...
         cgiters,errors2.cost,...
         cgiters,errors2.pi);
grid on
legend('S-viol.','S-cost', 'S-plan','N-viol.','N-cost','N-plan');
xlabel('errors over iterations (CG for Newton)');

figure(2);  
semilogy(time,errors.mu,  '--',...
         time,errors.cost,'--',...
         time,errors.pi,  '--'); hold on;
ax = gca; ax.ColorOrderIndex = 1;
semilogy(time2,errors2.mu,...
         time2,errors2.cost,...
         time2,errors2.pi);
grid on
legend('S-viol.','S-cost', 'S-plan','N-viol.','N-cost','N-plan');
xlabel('errors over run time [s]');
