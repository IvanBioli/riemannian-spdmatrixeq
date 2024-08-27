% Description: Test script for cost, gradient, Hessian and preconditioners
% in the metric Frobenius metric

clear all; clc;
%% GENERATE THE PROBLEM
% Note: the randomly generated problem not definite. However, this is not a
% problem for the purpose of this numerical test.
m = 1200; n = 1000; r_C = 10; fixed_rank = 10;
ell_half = 5; ell = 2 * ell_half;
for i = 1:ell_half
    A{i} = rand(m, m); A{i+ell_half} = A{i}';
    B{i} = rand(n, n); B{i+ell_half} = B{i}';
end
C.L = rand(m, r_C); C.R = rand(n, r_C);

%% DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
% Manifold
manifold = fixedrankembeddedfactory(m, n, fixed_rank);
problem.M = manifold;
problem.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem.grad = @(X, store) sylv_grad(X, A, B, C, store);
problem.hess = @(X, eta, store) sylv_hess(X, eta, A, B, C, store);

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
figure()
checkgradient(problem)

% Numerically check Hessian consistency
figure()
checkhessian(problem)
%% CHECK THE APPROXIMATE PRECONDITIONERS
m = 10; n = m; r = 2;
M = fixedrankembeddedfactory(m, n, r);
A = rand(m); A = A'*A + eye(m); B = rand(n); B = B*B' + eye(n);
X = M.rand();
eta = M.randvec(X);

% Exact solution via vectorization
A_{1} = A; A_{2} = speye(size(A)); B_{2} = B; B_{1} = speye(size(B));
C.L = rand(m, 1); C.R = rand(n, 1);
[g_exact, ~] = sylv_precon_kron(X, eta, A_, B_, C, store) ;
err_summary(g_exact, eta, X, A, B, M, "Kron")

% Fast solution
g = sylv_precon_AXDEXB(X, eta, A, B, struct());
err_summary(g, eta, X, A, B, M, "Approx Newton fast")

% Approximate solution via ADI
eA = eig(A); a = min(eA); b = max(eA);
eB = eig(B); c = min(eB); d = max(eB);
[p, q] = zolotarev_poles(100, a, b, c, d);
options.p = p; options.q = q;
g = sylv_precon_tangADI(X, eta, A, B, struct(), options);
err_summary(g, eta, X, A, B, M, "Approx Newton ADI")

%% AUXILIARY FUNCIONS
function M = convert(X)
% convert - Converts X to a full matrix

M = X.U * X.S * X.V';
end
function err = error_solve(g, eta, X, A, B, M)
% error_solve - Computes the relative residual norm

eta_ambient = M.tangent2ambient(X, eta);
g_ambient = convert(M.tangent2ambient(X, g));
err = norm(convert(M.tangent2ambient(X, M.proj(X, A * g_ambient + g_ambient * B))) - convert(eta_ambient), 'fro') / norm(convert(eta_ambient), 'fro');
end

function err_summary(g, eta, X, A, B, M, name)
% err_summary - Prints an error summary

U = X.U; V = X.V;
fprintf(name + "\n");
fprintf("\t err. \t %e\n", error_solve(g, eta, X, A, B, M));
res_U = (A - U * U' * A) * g.Up + g.Up * (V'*B*V) - eta.Up + (A - U * U' * A) * U * g.M;
fprintf("\t eqU: \t %e\n", norm(res_U, 'fro'))
res_V = (B - V * V' * B) * g.Vp + g.Vp * (U'*A*U) - eta.Vp + (B - V * V' * B) * V * g.M';
fprintf("\t eqV: \t %e\n", norm(res_V, 'fro'))
res_M = (U' * A * U) * g.M + g.M * (V' * B * V) - eta.M + (U' * A) * g.Up + ((V' * B) * g.Vp)';
fprintf("\t eqM: \t %e\n", norm(res_M, 'fro'))
end