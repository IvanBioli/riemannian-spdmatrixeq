function [g, store] = sylv_richardson_fADI(X, A, B, C, store, options)
% sylv_richardson_fADI -  Applies the Riemannian truncated
% preconditioned Richardson iteration with generalized-Sylvester-like
% preconditioner. The inverse of the generalized Sylvester operator is
% approximated by steps of fADI algorithm
%
% Syntax:
%   [g, store] = sylv_richardson_fADI(X, A, B, C, store, options)
%
% Inputs:
%   - X      : Current iterate
%   - A, B   : Cell arrays containing the coefficient matrices of the
%              multiterm linear matrix equation
%   - C      : Structure of the right-hand-side. It is a structure with
%              fields L and R, and the RHS is C.L * C.R'
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure. It is a structure with fields:
%                - fadi: structure with fields A, B, D, E, p, q. The matrices
%                  A, B, D, E define the preconditioning generalized Sylvester
%                  operator P(X) = A * X * B + E * X * D. The vectors p and
%                  q contain the ADI shift parameters.
%                  If D, E are not specified they are assumed to be the
%                  identity matrix
%                - trunc_opts: structure with truncation options for
%                  truncating the Euclidean gradient before applying fADI.
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
if isfield(options, "trunc_opts")
    g = truncation(g, options.trunc_opts);
end

% Apply fADI iteration
if ~isfield(options.fadi, 'D')
    g_precond = fadi_precond(g, options.fadi.A, options.fadi.B, ...
        options.fadi.p, options.fadi.q);
else
    g_precond = fadi_precond(g, options.fadi.A, options.fadi.B, ...
        options.fadi.p, options.fadi.q, options.fadi.D, options.fadi.E);
end
% Project to the tangent
g = project(X, g_precond);

end

function g = project(X, g)
% project - Projects the preconditioned Euclidean gradient to the tangent
% space

% Project to the tangent
ZV = g.U * (g.V' * X.V);
ZtU = g.V * (g.U' * X.U);
UtZV = X.U'*ZV;
g.M = UtZV;
g.Up = ZV  - X.U*UtZV;
g.Vp = ZtU - X.V*UtZV';
end
