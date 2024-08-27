function [g, store] = sylv_precon_kronnaive(X, H, A, B, C, store)
% sylv_precon_kronnaive - Applies the Riemannian preconditioner
% obtained inverting the projected Euclidean Hessian, i.e. solves the
% equation Proj_X(calA(g)) = H for g in Tang_X(M_r). The projected equation
% is solved by vectorizing the equations on the factors M, Up and Vp.
%
% Syntax:
%   [g, store] = sylv_precon_kronnaive(X, H, A, B, C, store)
%
% Inputs:
%   - X    : Current iterate
%   - H    : Tangent vector to which the preconditioner is applied
%   - A, B : Cell arrays containing the coefficient matrices of the
%            multiterm linear matrix equation
%   - C    : Structure of the right-hand-side. It is a structure with
%            fields L and R, and the RHS is C.L * C.R'
%   - store: Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - g    : The result of the application of the preconditioner to eta
%   - store: Updated store struct
%

% Precomputations
store = sylv_prepare_cost(X, A, B, C, store);
[n, r] = size(X.V);
In = speye(n); Ir = speye(r); Pi_nr = TvecMat(n, r);
m = length(A);
U = X.U; V = X.V;

% Assemble the matrix LHS
L = zeros(n^2, r^2 + 2 * n * r);
for i = 1:m
    L = L + [kron(store.BV(:, :, i), store.AU(:, :, i)), ...
        kron(store.BV(:, :, i), A{i}),...
        kron(B{i}, store.AU(:, :, i)) * Pi_nr];
end

% Assemble the (giant) system
LHS = [...
    kron(V', U') * L;...
    kron(V', In - U*U') * L; ...
    kron(U', In - V*V') * TvecMat(n, n) * L; ...
    zeros(r^2), kron(Ir, U'), zeros(r^2, n*r); ...
    zeros(r^2, r^2 + n*r), kron(Ir, V')];
RHS = [H.M(:); H.Up(:); H.Vp(:); zeros(r^2, 1); zeros(r^2, 1)];

% Solve and reshape back
hv = LHS \ RHS; %norm(LHS * hv - RHS, 'fro') / norm(RHS, 'fro')
g.M = reshape(hv(1:r^2), [r, r]);
g.Up = reshape(hv(r^2+1:r^2+n*r), [n, r]);
g.Vp = reshape(hv(r^2+n*r+1:end), [n, r]);

% Ensure that g.M and g.Vp verify the tangent constraints
g_new.M = g.M;
g_new.Vp = g.Vp - V*(V'*g.Vp);
g_new.Up = g.Up - U*(U'*g.Up);

% Throw warning to detect problem
tol = 1e-10;
err_Vp = norm(g.Vp - g_new.Vp, 'fro')/ norm(g_new.Vp, 'fro');
if err_Vp > tol
    warning("Rel. error on tangent condition for g.Vp:%e\n", err_Vp)
end
err_Up = norm(g.Up - g_new.Up, 'fro')/ norm(g_new.Up, 'fro');
if err_Up > tol
    warning("Rel. error on tangent condition for g.Up:%e\n", err_Up)
end

g = g_new;
end
