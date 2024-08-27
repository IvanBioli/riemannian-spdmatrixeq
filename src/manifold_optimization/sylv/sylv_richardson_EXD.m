function [g, store] = sylv_richardson_EXD(X, A, B, C, Amean_d, Gmean_d, store)
% sylv_richardson_EXD - Applies the Riemannian
% truncated preconditioned Richardson iteration with preconditioner
% P(X) = Amean * X  * Gmean.
%
% Syntax:
%   [g, store] = sylv_richardson_EXD(X, A, B, C, Amean_d, Gmean_d, store)
%
% Inputs:
%   - X      : Current iterate
%   - A, B   : Cell arrays containing the coefficient matrices of the
%              multiterm linear matrix equation
%   - C      : Structure of the right-hand-side. It is a structure with
%              fields L and R, and the RHS is C.L * C.R'
%   - Amean_d: Decomposition of the matrix Amean
%   - Gmean_d: Decomposition of the matrix Gmean
%   - store  : Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - g    : Preconditioned search direction
%   - store: Updated store struct
%

% Precomputations
store = sylv_prepare_cost(X, A, B, C, store);

% Get problem size
[m, r] = size(X.U); n = size(X.V, 1);
ell = length(A);

% Assemble L(X) - C in factored form
g = struct();
g.U = [reshape(bsxfun(@times, store.AU, diag(X.S)'), [m, r*ell]), -C.L];
g.V = [reshape(store.BV, [n, r*ell]), C.R];

% Solve using Amean and Gmean
g.U = Amean_d \ g.U;
g.V = Gmean_d \ g.V;

% Project to the tangent
g = project(X, g);

end

function g = project(X, g)
% project - Projects g to the tangent space

% Project to the tangent
ZV = g.U * (g.V' * X.V);
ZtU = g.V * (g.U' * X.U);
UtZV = X.U'*ZV;
g.M = UtZV;
g.Up = ZV  - X.U*UtZV;
g.Vp = ZtU - X.V*UtZV';
end
