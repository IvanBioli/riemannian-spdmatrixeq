function [g, store] = sylv_precon_EX(X, H, Amean, Amean_d, store)
% sylv_precon_EX - Applies to H the Riemannian
% preconditioner derived form the operator P(X) = Amean * X, i.e. solves
% the equation Proj_X(Amean * g) = H, with g in Tang_X(M_r).
%
% Syntax:
%   [g, store] = sylv_precon_EX(X, H, Amean, Amean_d, store)
%
% Inputs:
%   - X      : Current iterate
%   - H      : Tangent vector to which the preconditioner is applied
%   - Amean  : Left matrix in the preconditioner
%   - Amean_d: Decomposition of Amean
%   - store  : Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - g    : The result of the application of the preconditioner to H
%   - store: Updated store struct
%

% Precomputations
U = X.U;
AmeanU = Amean * U;
UtAmeanU = U' * AmeanU; UtAmeanU = 0.5 * (UtAmeanU + UtAmeanU');% Ensure that it is numerically symmetric

% Compute g.Vp
g.Vp = H.Vp / UtAmeanU;

% Compute g.M and g.Up
aux = Amean_d \ (H.Up + U * H.M);
g.M = U' * aux;
g.Up = aux - U * g.M;
end
