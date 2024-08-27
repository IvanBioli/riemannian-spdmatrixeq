function [h, store] = lyap_psd_hess_MXM(X, eta, A, M, N, B, store, options, Mass_d)
% lyap_psd_hess_MXM - Computes the Hessian-vector product
% between the Riemannian or projected Euclidean Hessian, in the metric
% induced by P(X) = Mass * X * Mass, and the tangent vector eta.
%
% Syntax:
%   [h, store] = lyap_psd_hess_MXM(X, eta, A, M, N, B, store, options, Mass_d)
%
% Inputs:
%   - X      : Current iterate
%   - eta    : Tangent vector at X
%   - A, M   : Coefficient matrices of the multiterm linear matrix equation
%   - N      : Cell array containing the other coefficient matrices of the
%              multiterm linear matrix equation
%   - B      : Factor of the righ-hand-side B*B'
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure. See mergeDefaultOptions.
%   - Mass_d : Decomposition of the matrix Mass
%
% Outputs:
%   - h    : Hessian-vector product
%   - store: Updated store structure.
%

if nargin < 8
    options = struct();
end
options = mergeDefaultOptions(options);

% Precomputations
V = X.V; d = diag(X.D);
[n, r] = size(V);
ell = length(N);
store = lyap_psd_prepare_cost(X, A, M, N, B, store);
is_M = ~isempty(M);

% Compute V' * A * V_eta, V' * M * V_eta and V' * N_i * V_eta
VtAVeta = store.AV' * eta.Vp;
if is_M
    VtMVeta = store.MV' * eta.Vp;
end
VtNVeta = pagemtimes(store.NV, "transpose", eta.Vp, "none");

% Compute Z * V_eta
if ~options.projehess
    ZVeta = (store.MV .* d') * VtAVeta + ...
        - sum(pagemtimes(bsxfun(@times, store.NV, d'), VtNVeta), 3) ...
        - B * (B' * eta.Vp);
    if is_M
        ZVeta = ZVeta + (store.AV .* d') * VtMVeta;
    end
end

% Compute \dot{Z} * V
NVMNVeta = zeros(n, r, ell);
etaV = V * eta.M + eta.Vp;
for i = 1:ell
    NVMNVeta(:, :, i) = N{i} * etaV;
end
ZdotV = pagemtimes(NVMNVeta, store.VtNV) + ...
    pagemtimes(store.NV, "none", VtNVeta, "transpose");
ZdotV = -sum(ZdotV, 3) + A * etaV * store.VtMV + store.MV * VtAVeta';
if is_M
    ZdotV = ZdotV + M * etaV * store.VtAV + store.AV * VtMVeta';
else
    ZdotV = ZdotV + etaV * store.VtAV;
end

% Compute the components of the Hessian
h.M = V' * ZdotV; h.M = symmetrize(h.M);
if ~options.projehess
    h.MVp = ZdotV + ZVeta .* (1 ./ d');
    h.MVp = h.MVp - X.MV * (V' * h.MVp);
else
    h.MVp = ZdotV - X.MV * h.M;
end
h.Vp = Mass_d \ h.MVp;
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