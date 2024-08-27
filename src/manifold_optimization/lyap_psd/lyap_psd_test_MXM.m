% Description: Test script for cost, gradient, Hessian and preconditioners
% in the metric induced by P(X) = Mass * X * Mass

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
Mass = rand(n, n); Mass = Mass * Mass'; Mass_d = decomposition(Mass);
% Manifold
manifold = sympsdfixedrank_metricMXM_factory(Mass, fixed_rank);
problem.M = manifold;
% Cost
problem.cost = @(X, store) lyap_psd_cost(X, A, M, N, B, store);
% Gradient
problem.grad = @(X, store) lyap_psd_grad_MXM(X, A, M, N, B, store, Mass_d);
% Hessian
problem.hess = @(X, eta, store) lyap_psd_hess_MXM(X, eta, A, M, N, B, store, struct(), Mass_d);
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
figure();
checkgradient(problem)

figure();
checkhessian(problem);
%% CHECK THE APPROXIMATE NEWTON PRECONDITIONER
n = 10; r = 2;
A = rand(n); A = A'*A; Mass = rand(n); Mass = Mass * Mass';

M = sympsdfixedrank_metricMXM_factory(Mass, r);
X = M.rand();
eta = M.randvec(X); eta.MVp = Mass * eta.Vp;
[g, store] = lyap_psd_precon_AXMMXA(X, eta, A, struct(), struct(), Mass);
fprintf("Approx. Newton fast:\n\terror: \t %e\n", error_solve_metricMXM(g, eta, X, A, Mass, M))

M = sympsdfixedrankembeddedfactory(n, r);
X = M.rand();
eta = M.randvec(X);
eAM = eig(A, Mass); a = min(eAM); b = max(eAM);
p = zolotarev_poles(1000, a, b);
options.p = p;
[g, store] = lyap_psd_precon_tangADI(X, eta, A, struct(), options, Mass);
fprintf("Approx. Newton ADI:\n\terror: \t %e\n", error_solve(g, eta, X, A, Mass, M))

function M = convert(X)
% convert - Converts X to a full matrix

M = X.V * X.D * X.V';
end

function err = error_solve(g, eta, X, A, Mass, M)
% error_solve - Computes the relative residual norm in the Frobenius norm

eta_ambient = M.tangent2ambient(X, eta);
g_ambient = convert(M.tangent2ambient(X, g));
err = norm(convert(M.tangent2ambient(X, M.proj(X,  A * g_ambient * Mass + Mass * g_ambient * A))) - convert(eta_ambient), 'fro') / norm(convert(eta_ambient), 'fro');
end
function err = error_solve_metricMXM(g, eta, X, A, Mass, M)
% error_solve_metricMXM - Computes the relative residual norm in the metric induced
% by P(X) = MXM

eta_ambient = M.tangent2ambient(X, eta);
g_ambient = convert(M.tangent2ambient(X, g));
err = norm(convert(M.tangent2ambient(X, M.proj(X,  Mass \ (A * g_ambient) + (g_ambient * A) / Mass))) - convert(eta_ambient)) / norm(convert(eta_ambient), 'fro');
end