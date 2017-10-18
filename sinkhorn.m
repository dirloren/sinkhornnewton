function [pi,v,w,errors,time] = sinkhorn(mu,nu,cost,gamma,maxiter,check,tol,varargin)
%[pi,v,w,errors,time,cgiters] = sinkhornn(mu,nu,cost,gamma,maxiter,check,tol,varargin)
%SINKHORN - calculate entropically regularized transport distance with Sinkhorn-Knopp method
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
%    v,w - Lagrange multipliers
%    error - struct containing various errors
%    time - time stamps of the iterations

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

K = exp(-cost/gamma); % exp(-cost/gamma)

% Initialization 
v = ones(N1,1);
w = ones(N2,1);
Kw = K*w;

errmu = zeros(maxiter,1);
time = zeros(maxiter,1);
if exist('piopt','var')
    errpi = zeros(maxiter,1);
end
if exist('optcost','var')
    errcost = zeros(maxiter,1);
end

tic % ready, set, go
for k = 1:maxiter
    
    % get errors and timing
    errmu(k) = norm(v.*(Kw)-mu,Inf);
    time(k) = toc;
    
    % check if something went wrong
    if isnan(errmu(k)) || isinf(errmu(k))
        error('Sinkhorn failed, some quantity is not finite anymore')
    end
    
    % assemble pkan if needed
    if exist('piopt','var') || exist('optcost','var')
        pi = v.*(K.*w');
    end
    if exist('piopt','var')
        errpi(k) = norm(pi(:)-piopt(:),1);
    end
    if exist('optcost','var')
        tcost = sum(sum(cost.*pi));
        errcost(k) = abs(tcost-optcost);
    end
    
    % check for output, convergence or termination
    if mod(k,check)==0 || max(errmu(k))<tol || k==maxiter
        if ~(exist('piopt','var') || exist('optcost','var'))
            pi = v.*(K.*w');
        end
        
        % get errors in last step
        if ~exist('optcost','var')
            tcost = sum(sum(cost.*pi));
        end
        
        % output
        negent = sum(sum(pi.*log(pi+eps)));
        fprintf('iter: %d, mismatch in mu: %2.2e, tcost: %2.2e, negent: %2.2e, obj: %2.2e\n',...
            k,errmu(k),tcost,negent,tcost+gamma*negent)
        
        % if converged, exit here
        if  max(errmu(k))<tol
            errors.mu = errmu(1:k);
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
            error('Sinkhorn did not converge within the maximum number of iterations')
        end
    end
    % Sinkhorn-Knopp
    % iteratively scale rows and columns
    v = mu./(Kw);
    w = nu./(K'*v);
    Kw = K*w;    
        
end
toc
