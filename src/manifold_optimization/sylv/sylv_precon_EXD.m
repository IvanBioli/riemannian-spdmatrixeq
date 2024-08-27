function [g, store] = sylv_precon_EXD(X, H, A, A_d, B, B_d, store)
% sylv_precon_EXD - Applies to H the
% Riemannian preconditioner derived form the operator P(X) = AXB, i.e.
% solves the equation Proj_X(A * g * B) = H, with g in Tang_X(M_r).
%
% Syntax:
%   [g, store] = sylv_precon_EXD(X, H, A, A_d, B, B_d, store)
%
% Inputs:
%   - X    : Current iterate
%   - H    : Tangent vector to which the preconditioner is applied
%   - A    : Left matrix in the preconditioner
%   - A_d  : Decomposition of A
%   - B    : Right matrix in the preconditioner
%   - B_d  : Decomposition of B
%   - store: Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - g    : The result of the application of the preconditioner to H
%   - store: Updated store struct
%


% Precomputations
U = X.U; V = X.V;

AU = A * U;
UtAU = U' * AU;
UtAU = 0.5 * (UtAU + UtAU'); % Ensure that it is numerically symmetric
UtAU_d = decomposition(UtAU, 'chol');

BV = B * V;
VtBV = V' * BV;
VtBV = 0.5 * (VtBV + VtBV'); % Ensure that it is numerically symmetric
VtBV_d = decomposition(VtBV, 'chol');

% Compute g.Up
aux_U = A_d \ (H.Up + U * H.M) / VtBV_d;
g.Up = aux_U - U * (U' * aux_U);

% Compute g.Vp
aux_V = B_d \ (H.Vp + V * H.M') / UtAU_d;
g.Vp = aux_V - V * (V' * aux_V);

% Compute g.M
g.M = UtAU_d \ (H.M - AU' * g.Up * VtBV - UtAU * g.Vp' * BV) / VtBV_d;
end
