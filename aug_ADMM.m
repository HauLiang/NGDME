function [x, history] = aug_ADMM(A, D, b, lambda, rho, alpha)
% Solves general lasso problem via aug-ADMM
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda || z ||_1
%   subject to z = Dx
%
%   The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
% More information can be found in these papers linked at:
% https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1114491
% https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
%
% Author: Hao Liang
% Last modified by: 21/04/13
%

% Global constants and defaults
MAX_ITER = 1000;
ABSTOL   = 1e-4;
RELTOL   = 1e-4;

% Data preprocessing
[~, n] = size(A);
[m1, ~] = size(D);

% save a matrix-vector multiply
Atb = A'*b;
AtA = A'*A;
DtD = D.'*D;

% aug-ADMM solver
x = zeros(n,1);
z = zeros(m1,1);
u = zeros(m1,1);
u_pre = zeros(m1,1);

% set the additional variables
T = 1*(D'*D);

for k = 1:MAX_ITER

    % x-update
    q = Atb - D'*(2*u-u_pre)+rho/2*(T+T')*x;    % temporary value
    x = (AtA + rho*(DtD)) \ q;

    % z-update using soft threshold operator
    zold = z;
    x_hat = alpha*D*x + (1 - alpha)*zold;
    z = shrinkage(x_hat + u, lambda/rho);
    
    % u-update
    u_pre = u;
    u = u + rho*(x_hat-z);

    % diagnostics, reporting, termination checks
    history.objval(k)  = abs(objective(A, b, lambda, x, z));

    history.r_norm(k)  = abs(norm(D*x - z));
    history.s_norm(k)  = abs(norm(-rho*(z - zold)));

    history.eps_pri(k) = abs(sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z)));
    history.eps_dual(k)= abs(sqrt(n)*ABSTOL + RELTOL*norm(rho*u));

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

end

end

% objective function
function p = objective(A, b, lambda, x, z)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );
end

% soft-threshold operator
function z = shrinkage(x, kappa)
    temp = abs(x);
    if temp < kappa
        z = zeros(size(x));
    else
        z = (temp>kappa).* (temp-kappa).*(x./(temp+eps));
    end
end

