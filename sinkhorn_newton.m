function [pi,alpha,beta,errors,time,cgiters] = sinkhorn_newton(mu,nu,cost,gamma,maxiter,check,tol,varargin)
%[pi,alpha,beta,errors,time,cgiters] = sinkhorn_newton(mu,nu,cost,gamma,maxiter,check,tol,varargin)
%SINKHORN_NEWTON - calculate entropically regularized transport distance with Sinkhorn-Newton method
%
% Inputs:
%    mu, nu - non-negative vectors with equal sum and length N1 and N2, resp.
%    cost - cost matrix of size N1 x N2
%    gamma - positive regularization parameter
%    maxiter - maximum number of iterations
%    check - output frequency in sinkhorn_newton
%    tol - tolerance for optimality check in sinkhorn_newton
% Optional input:
%  'optplan' -  optimal plan (sime size as pi), used to keep track of error in plan
%  'optcost' - optimal transport cost sum(sum(c.*piopt)), used to keep track
%     of error in cost
% Outputs:
%    pi - optimal transport plan
%    alpha, beta - basically log of Lagrange multipliers
%    error - struct containing various errors
%    time - time stamps of the iterations
%    cgiters - number of cg iterations performed up the respective Newton
%       step

% Author: Dirk Lorenz
% email: d.lorenz@tu-bs.de
% Website: https://www.tu-braunschweig.de/iaa/personal/lorenz
% October 2017; Last revision: 17-October-2017

% parse varargin
for i=1:2:(length(varargin)-1)
    switch varargin{i}
        case 'optplan'
            piopt = varargin{i+1};
        case 'optcost'
            optcost = varargin{i+1};
    end
end

N1 = length(mu);
N2 = length(nu);

% Initialization 
alpha = zeros(N1,1);
beta = zeros(N2,1);

errnu = zeros(maxiter,1);
errmu = zeros(maxiter,1);
cgiters = zeros(maxiter,1);
time = zeros(maxiter,1);
if exist('piopt','var')
    errpi = zeros(maxiter,1);
end
if exist('optcost','var')
    errcost = zeros(maxiter,1);
end
alliters = 0;
e1 = ones(N1,1);
e2 = ones(N2,1);
dx = zeros(N1+N2,1);

tic; % ready, set, go
for k = 1:maxiter+1
    
    % assemble plan and compute marginals
    pi = exp((-cost - alpha - beta')/gamma);
    pie = sum(pi,2);
    pite = sum(pi,1)';
    
    % get errors and timing
    errmu(k) = norm(pie-mu,Inf);
    errnu(k) = norm(pite-nu,Inf);
    time(k) = toc;
    if exist('piopt','var')
        errpi(k) = norm(pi(:)-piopt(:),1);
    end
    if exist('optcost','var')
        tcost = sum(sum(cost.*pi));
        errcost(k) = abs(tcost-optcost);
    end
    
    % check if something went wrong
    if isnan(errmu(k)) || isinf(errmu(k))
        error('sinkhorn_newton failed, some quantity is not finite anymore')
    end
    
    % check for output, convergence or termination
    if mod(k,check)==0 || max(errnu(k),errmu(k))<tol || k==maxiter
        % get errors in last step
        if ~exist('optcost','var')
            tcost = sum(sum(cost.*pi));
        end
        
        % output
        negent = sum(sum(pi.*log(pi+eps)));
        fprintf('iter: %d, mismatch in mu: %2.2e, mismatch in nu: %2.2e, tcost: %2.2e, negent: %2.2e, obj: %2.2e\n',...
            k,errmu(k),errnu(k),tcost,negent,tcost+gamma*negent)
        
        % if converged, exit here
        if max(errnu(k),errmu(k))<tol
            errors.mu = errmu(1:k);
            errors.nu = errnu(1:k);
            cgiters = cgiters(1:k);
            time = time(1:k);
            if exist('piopt','var')
                errors.pi = errpi(1:k);
            end
            if exist('optcost','var')
                errors.cost = errcost(1:k);
            end
            toc
            return
        end
        
        % if maxiter reched, exit here
        if k==maxiter
            error('sinkhorn_newton did not converge within the maximum number of iterations')
        end
    end
    % Newton step
    H = @(x) [pie.*x(1:N1) + pi*x(1+N1:end); pi'*x(1:N1) + pite.*x(1+N1:end)];
    g = gamma*[mu-pie;nu-pite];
    P = @(x) x./([pie;pite]+1e-12);
    [dx,flag,relres,iters,res] = pcg(H,g,tol,max(round((N1+N2)/24),10),P,[],dx);
    % reshift for equal sums in alpha beta
    dotProd = (sum(dx(1:N1))-sum(dx(N1+1:end)))/(N1+N2);
    dx(1:N1) = dx(1:N1) - dotProd; dx(N1+1:end) = dx(N1+1:end) + dotProd;
    alliters = alliters+iters;
    cgiters(k+1) = alliters;
    alpha = alpha - dx(1:N1);
    beta = beta - dx(N1+1:end);
end
toc
