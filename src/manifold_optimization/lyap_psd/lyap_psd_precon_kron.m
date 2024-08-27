function [g, store] = lyap_psd_precon_kron(X, H, A, M, N, store)
% lyap_psd_precon_kron - Applies the Riemannian preconditioner
% obtained inverting the projected Euclidean Hessian, i.e. solves the
% equation Proj_X(calA(g)) = H for g in Tang_X(Snr++). The projected equation
% is solved by vectorizing the equations on the factors M and Vp.
%
% Syntax:
%   [g, store] = lyap_psd_precon_kron(X, H, A, M, N, store)
%
% Inputs:
%   - X    : Current iterate
%   - H    : Tangent vector at X
%   - A    : Coefficient matrix of the multiterm linear matrix equation
%   - M    : Coefficient matrix of the multiterm linear matrix equation.
%            Neglected, i.e. assumed to be the identity
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

AV = A*V;
VtA = V' * A;
VtAV = VtA * V;
VtN = cell(m, 1); VtNV = cell(m, 1); NV = cell(m, 1);
for i = 1:m
    VtN{i} = V' * N{i};
    VtNV{i} = VtN{i} * V;
    NV{i} = N{i} * V;
end

% Assemble the matrix L

LHS_hh = kron(Ir, VtAV) + kron(VtAV, Ir);
LHS_hv = kron(Ir, VtA) + kron(VtA, Ir) * Pi;
LHS_vh = kron(Ir, AV - V * VtAV);
LHS_vv = kron(Ir, A - V*VtA) + kron(VtAV, In);

for i = 1:m
    LHS_hh = LHS_hh + kron(VtNV{i}, VtNV{i});
    LHS_hv = LHS_hv + kron(VtNV{i}, VtN{i}) + kron(VtN{i}, VtNV{i}) * Pi;
    LHS_vh = LHS_vh + kron(VtNV{i}, NV{i} - V*VtNV{i});
    LHS_vv = LHS_vv + kron(VtNV{i}, N{i} - V*VtN{i}) + kron(VtN{i}, NV{i} - V*VtNV{i}) * Pi;
end

% Assemble the (giant) system
LHS = [...
    LHS_hh, LHS_hv;...
    LHS_vh, LHS_vv...
    ];
RHS = [H.M(:); H.Vp(:)];

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
