function [h, store] = sylv_hess_EXD(X, eta, A, B, C, store, E_d, D_d, options)
% sylv_hess_EXD - Computes the Hessian-vector product
% between the Riemannian or projected Euclidean Hessian, in the metric
% induced by P(X) = EXD, and the tangent vector eta.
%
% Syntax:
%   [h, store] = sylv_hess_EXD(X, eta, A, B, C, store, E_d, D_d, options)
%
% Inputs:
%   - X      : Current iterate
%   - eta    : Tangent vector at X
%   - A, B   : Cell arrays containing the coefficient matrices of the
%              multiterm linear matrix equation
%   - C      : Structure of the right-hand-side. It is a structure with
%              fields L and R, and the RHS is C.L * C.R'
%   - store  : Manopt's StoreDB struct for the current iterate
%   - E_d    : Decomposition of the matrix E
%   - D_d    : Decomposition of the matrix D
%   - options: Options structure. See mergeDefaultOptions.
%
% Outputs:
%   - h    : Hessian-vector product
%   - store: Updated store structure.
%


if nargin < 9
    options = struct();
end
options = mergeDefaultOptions(options);

% Precomputations
U = X.U; V = X.V;
[m, r] = size(U); [n, ~] = size(V);
ell = length(A); s = diag(X.S);
store = sylv_prepare_cost(X, A, B, C, store);

% Compute B_i * V_eta and A_i * U_eta
AUeta = zeros(m, r, ell); AUMeta  = zeros(m, r, ell); UMeta  = U * eta.M;
BVeta = zeros(n, r, ell); BVMetat = zeros(n, r, ell); VMetat = V * eta.M';
for i = 1:ell
    AUeta(:, :, i) = A{i} * eta.Up; AUMeta(:, :, i)  = A{i} * UMeta;
    BVeta(:, :, i) = B{i} * eta.Vp; BVMetat(:, :, i) = B{i} * VMetat;
end

% Compute Z * V_eta and Z' * U_eta
if ~options.projehess
    ZVeta  = sum(pagemtimes(bsxfun(@times, store.AU, s'), ...
        pagemtimes(store.BV, "transpose", eta.Vp, "none")), 3) ...
        - C.L * (C.R' * eta.Vp);
    ZtUeta = sum(pagemtimes(bsxfun(@times, store.BV, s'), ...
        pagemtimes(store.AU, "transpose", eta.Up, "none")), 3) ...
        - C.R * (C.L' * eta.Up);
end

% Compute \dot{Z} * V and \dot{Z}' * U
ZdotV  = pagemtimes(AUeta + AUMeta, "none", store.VtBV, "transpose") + ...
    pagemtimes(store.AU, pagemtimes(BVeta, "transpose", V, "none"));
ZdotV  = sum(ZdotV, 3);
ZdottU = pagemtimes(BVeta + BVMetat, "none", store.UtAU, "transpose") + ...
    pagemtimes(store.BV, pagemtimes(AUeta, "transpose", U, "none"));
ZdottU = sum(ZdottU, 3);

% Compute the components of the Hessian
h.M = U' * ZdotV;
if ~options.projehess
    s_reg = sqrt(s.^2 + eps^2);
    h.EUp = ZVeta .* (1 ./ s_reg') + ZdotV; h.EUp = h.EUp - X.EU * (U' * h.EUp);
    h.DVp = ZtUeta .* (1 ./ s_reg') + ZdottU; h.DVp = h.DVp - X.DV * (V' * h.DVp);
else
    h.EUp = ZdotV - X.EU * h.M;
    h.DVp = ZdottU - X.DV * h.M';
end
h.Up = E_d \ h.EUp;
h.Vp = D_d \ h.DVp;
end


function opts = mergeDefaultOptions(opts)
% mergeDefaultOptions - Merges the option structure given by the user with
%   the default option structure. User-defined options have precedence
%   over the default ones.
%
% Inputs:
%   - opts: Options structure with field opts.optionname, where optionname
%       is one of the following and the default value is indicated
%       between parentheses:
%           - projehess (false): whether to compute only the projected
%             Euclidean Hessian (i.e. do not compute the curvature term)


default_opts.projehess = false;
opts = mergeOptions(default_opts, opts);
end
