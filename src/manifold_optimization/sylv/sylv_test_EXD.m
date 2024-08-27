% Description: Test script for cost, gradient, Hessian and preconditioners
% in the metric induced by P(X) = E * X * D

clear all; clc;
%% GENERATE THE PROBLEM
% Note: the randomly generated problem is symmetric, but not definite (in
% general). However, this is not a problem for the purpose of this
% numerical test
n = 1000; m = 800; r_C = 5; fixed_rank = 10;
ell_half = 5; ell = 2 * ell_half;
for i = 1:ell_half
    A{i} = rand(m, m); A{i+ell_half} = A{i}';
    B{i} = rand(n, n); B{i+ell_half} = B{i}';
end
C.L = rand(m, r_C); C.R = rand(n, r_C);
E = rand(m, m); E = E*E';
D = rand(n, n); D = D*D';

%% DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
% Manifold
manifold = fixedrank_metricEXD_factory(E, D, fixed_rank);
problem.M = manifold;
problem.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem.grad = @(X, store) sylv_grad_EXD(X, A, B, C, store, decomposition(E), decomposition(D));
problem.hess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, decomposition(E), decomposition(D));

%% CHECK THE COST
X = problem.M.rand(); X_mat = manifold.triplet2matrix(X);
store = struct(); store.shared = struct();
[cost, store] = problem.cost(X, store);
fprintf("cost: %f\n", cost)
LX = zeros(m, n);
for i = 1:ell
    LX = LX + A{i} * X_mat * B{i}';
end
cost_exact = 0.5 * mat_inner(LX, X_mat) - mat_inner(X_mat, C.L * C.R');
fprintf("cost_exact: %f\n", cost_exact)

% Numerically check the cost
fprintf("rel err: %e\n", abs(cost - cost_exact) / abs(cost_exact))

%% CHECK THE GRADIENT
% Numerically check gradient consistency
figure();
checkgradient(problem)

% Numerically check Hessian consistency
figure();
checkhessian(problem);
%% CHECK THE APPROXIMATE NEWTON PRECONDITIONER
m = 10; n = 8; r = 2;
A = rand(m); A = A'*A; B = rand(n); B = B*B';
E = rand(m); E = E*E'; D = rand(n); D = D*D';

M = fixedrank_metricEXD_factory(E, D, r);
X = M.rand();
eta = M.randvec(X); eta.EUp = E * eta.Up; eta.DVp = D * eta.Vp;
[g, store] = sylv_precon_AXDEXB(X, eta, A, B, struct(), struct(), D, E);
fprintf("Error Newton fast: \t %e\n", error_solve_metricEXD(g, eta, X, A, B, D, E, M))

M = fixedrankembeddedfactory(m, n, r);
X = M.rand();
eta = M.randvec(X);
eAE = eig(A, E); a = min(eAE); b = max(eAE);
eBD = eig(B, D); c = min(eBD); d = max(eBD);
[p, q] = zolotarev_poles(1000, a, b, c, d);
options.p = p; options.q = q;
g = sylv_precon_tangADI(X, eta, A, B, struct(), options, E, D);
fprintf("Error ADI: \t %e\n", error_solve(g, eta, X, A, B, D, E, M));

function M = convert(X)
% convert - Converts X to a full metric

M = X.U * X.S * X.V';
end

function err = error_solve_metricEXD(g, eta, X, A, B, D, E, M)
% error_solve_metricEXD - Computes the relative residual norm in the metric induced
% by P(X) = EXD

eta_ambient = M.tangent2ambient(X, eta);
g_ambient = convert(M.tangent2ambient(X, g));
err = norm(convert(M.tangent2ambient(X, M.proj(X,  E \ (A * g_ambient) + (g_ambient * B) / D))) - convert(eta_ambient), 'fro') / norm(convert(eta_ambient), 'fro');
end

function err = error_solve(g, eta, X, A, B, D, E, M)
% error_solve - Computes the relative residual norm in the Frobenius norm

eta_ambient = M.tangent2ambient(X, eta);
g_ambient = convert(M.tangent2ambient(X, g));
err = norm(convert(M.tangent2ambient(X, M.proj(X,  A * g_ambient * D + E * g_ambient * B))) - convert(eta_ambient), 'fro') / norm(convert(eta_ambient), 'fro');
end
