function [g, store] = sylv_precon_kron(X, H, A, B, C, store)
% sylv_precon_kron - Applies the Riemannian preconditioner
% obtained inverting the projected Euclidean Hessian, i.e. solves the
% equation Proj_X(calA(g)) = H for g in Tang_X(M_r). The projected equation
% is solved by vectorizing the equations on the factors M, Up and Vp.
%
% Syntax:
%   [g, store] = sylv_precon_kron(X, H, A, B, C, store)
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
%   - g    : The result of the application of the preconditioner to H
%   - store: Updated store struct
%

% Precomputations
store = sylv_prepare_cost(X, A, B, C, store);
[m, r] = size(X.U); [n, ~] = size(X.V);
Ir = speye(r); Pi_nr = TvecMat(n, r); Pi_mr = TvecMat(m, r); Pi_rr = TvecMat(r,r);
ell = length(A);
U = X.U; V = X.V;
AU = cell(ell, 1); BV = cell(ell, 1); UtAU = cell(ell, 1); VtBV = cell(ell, 1);
for i = 1:ell
    AU{i} = A{i} * U;
    BV{i} = B{i} * V;
    UtAU{i} = U' * AU{i};
    VtBV{i} = V' * BV{i};
end

% Assemble the LHS blocks
LHS_mm = zeros(r^2); LHS_mu = zeros(r^2, m*r); LHS_mv = zeros(r^2, n*r);
LHS_um = zeros(m*r, r^2); LHS_uu = zeros(m*r); LHS_uv = zeros(m*r, n*r);
LHS_vm = zeros(n*r, r^2); LHS_vu = zeros(n*r, m*r); LHS_vv = zeros(n*r);
for i = 1:ell
    LHS_mm = LHS_mm + kron(VtBV{i}, UtAU{i});
    LHS_mu = LHS_mu + kron(VtBV{i}, AU{i}');
    LHS_mv = LHS_mv + kron(BV{i}', UtAU{i});
    
    LHS_um = LHS_um + kron(VtBV{i}, AU{i} - U*UtAU{i});
    LHS_uu = LHS_uu + kron(VtBV{i}, A{i} - U*AU{i}');
    LHS_uv = LHS_uv + kron(BV{i}', AU{i} - U*UtAU{i});
    
    LHS_vm = LHS_vm + kron(UtAU{i}, BV{i} - V*VtBV{i});
    LHS_vu = LHS_vu + kron(AU{i}', BV{i} - V*VtBV{i});
    LHS_vv = LHS_vv + kron(UtAU{i}, B{i} - V*BV{i}');
end

% Assemble the (giant) system
LHS = [...
    LHS_mm,             LHS_mu,             LHS_mv * Pi_nr;...
    LHS_um,             LHS_uu,             LHS_uv * Pi_nr; ...
    LHS_vm * Pi_rr,     LHS_vu * Pi_mr,     LHS_vv; ...
    zeros(r^2), kron(Ir, U'), zeros(r^2, n*r); ...
    zeros(r^2, r^2 + m*r), kron(Ir, V')];
RHS = [H.M(:); H.Up(:); H.Vp(:); zeros(r^2, 1); zeros(r^2, 1)];

% Solve and reshape back
hv = LHS \ RHS; %norm(LHS * hv - RHS, 'fro') / norm(RHS, 'fro')
g.M = reshape(hv(1:r^2), [r, r]);
g.Up = reshape(hv(r^2+1:r^2+m*r), [m, r]);
g.Vp = reshape(hv(r^2+m*r+1:end), [n, r]);

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
