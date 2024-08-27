function [g, store] = sylv_richardson_sylvesterexact(X, A, B, C, store, options)
% sylv_richardson_sylvesterexact - Applies the Riemannian truncated
% preconditioned Richardson iteration with generalized-Sylvester-like
% preconditioner. The preconditioner is applied exactly, i.e. the Sylvester
% equation is solved using exact and not iterative methods.
%
% Syntax:
%   [g, store] = sylv_richardson_sylvesterexact(X, A, B, C, store, options)
%
% Inputs:
%   - X      : Current iterate
%   - A, B   : Cell arrays containing the coefficient matrices of the
%              multiterm linear matrix equation
%   - C      : Structure of the right-hand-side. It is a structure with
%              fields L and R, and the RHS is C.L * C.R'
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure. It is a structure with fields A, B, D, E.
%              The (generalized) Sylvester operator used for preconditioning is
%              P(X) = options.A * X * options.D + options.E * X * options.B
%              If D, E are not specified they are assumed to be
%              the identity matrix
%
% Outputs:
%   - g    : Preconditioned search direction
%   - store: Updated store struct
%

% Precomputations
store = sylv_prepare_cost(X, A, B, C, store);
[m, r] = size(X.U); n = size(X.V, 1);
ell = length(A);

% Assemble L(X) - C in factored form
g = struct();
g.U = [reshape(bsxfun(@times, store.AU, diag(X.S)'), [m, r*ell]), -C.L];
g.V = [reshape(store.BV, [n, r*ell]), C.R];
g = g.U * g.V';

% Apply BartelStewart
if ~isfield(options, 'D')
    g_precond = sylvester(options.A, options.B, full(g));
else
    g_precond = bartelsStewart(options.A, options.D, options.E, options.B, full(g), false, false);
end
% Project to the tangent
g = project(X, g_precond);

end

function gproj = project(X, g)
% project - Projects the preconditioned Euclidean gradient to the tangent
% space

% Project to the tangent
ZV = g * X.V;
ZtU = g' * X.U;
UtZV = X.U'*ZV;
gproj.M = UtZV;
gproj.Up = ZV  - X.U*UtZV;
gproj.Vp = ZtU - X.V*UtZV';
end
