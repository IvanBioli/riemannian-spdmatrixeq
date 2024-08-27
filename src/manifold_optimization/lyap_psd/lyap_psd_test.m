% Description: Test script for cost, gradient, Hessian and preconditioners

clear all; clc;
%% GENERATE THE PROBLEM
n = 1000; r_B = 10; fixed_rank = 2;
m = 5;
A = rand(n, n); A = A*A';
M = rand(n, n); M = M*M';
N = cell(m, 1);
for i = 1:m
    N{i} = sqrt(1 / (2*m)) * rand(n, n); N{i} = N{i}*N{i}';
end
B = rand(n, r_B);

%% DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
% Manifold
manifold = sympsdfixedrankembeddedfactory(n, fixed_rank);
problem.M = manifold;
% Cost
problem.cost = @(X, store) lyap_psd_cost(X, A, M, N, B, store);
% Gradient
problem.grad = @(X, store) lyap_psd_grad(X, A, M, N, B, store);
% Hessian
problem.hess = @(X, eta, store) lyap_psd_hess(X, eta, A, M, N, B, store);
%% CHECK THE COST
X = problem.M.rand(); X_mat = manifold.couple2matrix(X);
store = struct(); store.shared = struct();
[cost, store] = problem.cost(X, store);
fprintf("cost: %f\n", cost)
LX = A * X_mat * M + M * X_mat * A;
for i = 1:m
    LX = LX - N{i} * X_mat * N{i};
end
cost_exact = 0.5 * mat_inner(LX, X_mat) - mat_inner(X_mat, B*B');
fprintf("cost_exact: %f\n", cost_exact)

% Numerically check the cost
fprintf("rel err: %e\n", abs(cost - cost_exact) / abs(cost_exact))

%% CHECK THE GRADIENT AND THE HESSIAN
figure()
checkgradient(problem)
figure()
checkhessian(problem)
%% CHECK THE APPROXIMATE NEWTON PRECONDITIONER
n = 10; r = 2;
manifold = sympsdfixedrankembeddedfactory(n, r);
A = rand(n); A = A'*A;
M = eye(n);
X = manifold.rand(); V = X.V;
eta = manifold.randvec(X);

N = {};
[g, ~] = lyap_psd_precon_kronnaive(X, eta, A, M, N, struct()) ;
err_summary(g, eta, X, A, M, manifold, "Naive Kron")

N = {};
[g, ~] = lyap_psd_precon_kron(X, eta, A, M, N, struct()) ;
err_summary(g, eta, X, A, M, manifold, "Kron")

[g, store] = lyap_psd_precon_AXMMXA(X, eta, A, struct());
err_summary(g, eta, X, A, M, manifold, "Approx Newton fast")

eA = eig(A); a = min(eA); b = max(eA);
options.p = zolotarev_poles(1000, a, b);
[g, store] = lyap_psd_precon_tangADI(X, eta, A, struct(), options);
err_summary(g, eta, X, A, M, manifold, "Approx Newton ADI")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = convert(X)
% convert - Converts X to a full metric

M = X.V * X.D * X.V';
end
function err = error_solve(g, eta, X, A, M, manifold)
% error_solve - Computes the relative residual norm in the Frobenius norm

eta_ambient = manifold.tangent2ambient(X, eta);
g_ambient = convert(manifold.tangent2ambient(X, g));
err = norm(convert(manifold.tangent2ambient(X, manifold.proj(X, A * g_ambient * M + M * g_ambient * A))) - convert(eta_ambient), 'fro') / norm(convert(eta_ambient), 'fro');
end
function err_summary(g, eta, X, A, M, manifold, name)
% err_summary - Prints an error summary

V = X.V;
fprintf(name + "\n");
fprintf("\t err. \t %e\n", error_solve(g, eta, X, A, M, manifold));
res_V = (A - V * V' * A) * g.Vp + g.Vp * (V'*A*V) - eta.Vp + (A - V * V' * A) * V * g.M;
fprintf("\teqV: \t %e\n", norm(res_V, 'fro'))
res_M = (V' * A * V) * g.M + g.M * (V' * A * V) - eta.M + (V' * A) * g.Vp + ((V' * A) * g.Vp)';
fprintf("\teqM: \t %e\n", norm(res_M, 'fro'))
end