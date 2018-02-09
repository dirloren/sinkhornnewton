offset = 0.01;
maxiter = 100000;
check = inf;
tol = 1e-12;

% MNIST images, can be obtained from http://yann.lecun.com/exdb/mnist/
images = loadMNISTImages('train-images.idx3-ubyte');  
p = [4378 37]; % index of images to compute transport between

mu = images(:, p(1)) + offset;
mu = mu / sum(mu);
nu = images(:, p(2)) + offset;
nu = nu / sum(nu);

cost = MNISTGroundMetric(length(mu), 1 / (sqrt(length(mu)) - 1)).^2;

close all; 
for gamma = [1 .1 .01 .005] * median(cost(:))
    [pi,errors,time,cgiters] = sinkhorn_newton_primal(mu,nu,cost,gamma,maxiter,check,tol);
    
    figure(1); semilogy(cgiters,errors.mu); hold on; grid on; 
    figure(2); semilogy(time,errors.mu); hold on; grid on; drawnow
end

figure(1)
legend('\gamma = q_{50}','\gamma = 0.1q_{50}','\gamma = 0.01q_{50}','\gamma = 0.005 q_{50}');
xlabel(['\epsilon = ' num2str(offset) ' (CG iterations)']);

figure(2)
legend('\gamma = q_{50}','\gamma = 0.1q_{50}','\gamma = 0.01q_{50}','\gamma = 0.005 q_{50}');
xlabel(['\epsilon = ' num2str(offset) ' (run time in seconds)']);
