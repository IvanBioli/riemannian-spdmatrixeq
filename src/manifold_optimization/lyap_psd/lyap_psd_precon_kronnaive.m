function [g, store] = lyap_psd_precon_kronnaive(X, H, A, M, N, store)
% lyap_psd_precon_kronnaive - Applies the Riemannian preconditioner
% obtained inverting the projected Euclidean Hessian, i.e. solves the
% equation Proj_X(calA(g)) = H for g in Tang_X(Snr++). The projected equation
% is solved by vectorizing the equations on the factors M and Vp.
%
% Syntax:
%   [g, store] = lyap_psd_precon_kronnaive(X, H, A, M, N, store)
%
% Inputs:
%   - X    : Current iterate
%   - H    : Tangent vector at X
%   - A, M : Coefficient matrices of the multiterm linear matrix equation
%   - N    : Cell array containing the other coefficient matrices of the
%            multiterm linear matrix equation
%   - store: Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - g    : The result of the application of the preconditioner to H
%   - store: Updated store struct
%



% Precomputations
[n, r] = size(X.V);
In = speye(n); Ir = speye(r); Pi = TvecMat(n, r);
m = length(N);
V = X.V;

% Assemble the matrix L
L = kron(A, M) + kron(M, A);
for i = 1:m
    L = L - kron(N{i}, N{i});
end
L = L * [kron(V,V), kron(V, In) + kron(In, V) * Pi];

% Assemble the (giant) system
LHS = [...
    kron(V', V') * L;...
    kron(V', In - V*V') * L; ...
    zeros(r^2), kron(Ir, V'); ...
    -speye(r^2) + TvecMat(r, r), zeros(r^2, n*r)];
RHS = [H.M(:); H.Vp(:); zeros(r^2, 1); zeros(r^2, 1)];

% Solve and reshape back
hv = LHS \ RHS;
g.M = reshape(hv(1:r^2), [r, r]);
g.Vp = reshape(hv(r^2+1:end), [n, r]);

% Ensure that g.M and g.Vp verify the tangent constraints
g_new.M = 0.5 * (g.M + g.M');
g_new.Vp = g.Vp - V*(V'*g.Vp);

% Throw warning to detect problem
tol = 1e-7;
err_M = norm(g.M - g_new.M, 'fro') / norm(g_new.M, 'fro');
err_Vp = norm(g.Vp - g_new.Vp, 'fro')/ norm(g_new.Vp, 'fro');
if err_Vp > tol
    warning("Rel. error on tangent condition for g.Vp:%e\n", err_Vp)
end
if err_M > tol
    warning("Rel. error on tangent condition for g.M:%e\n", err_M)
end

g = g_new;
end
