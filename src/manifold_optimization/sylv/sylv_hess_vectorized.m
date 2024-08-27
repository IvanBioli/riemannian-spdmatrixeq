function [Heta, store] = sylv_hess_vectorized(X, eta, A, B, C, store, M)
% sylv_hess_vectorized - Computes the Hessian-vector product
% between the Riemannian or projected Euclidean Hessian and the tangent
% vector eta, where everything is vectorized after passing in coordinates
% with respect to a basis of the tangent space. Not useful for practical
% instances, only for computing the eigenvalues of the Riemannian Hessian,
% projected Euclidean Hessian and curvature term.
%
% Syntax:
%   [Heta, store] = sylv_hess_vectorized(X, eta, A, B, C, store, M)

% Add basis and Hessians
if ~isfield(store, "basis")
    store.basis = tangentorthobasis(M, X);
end
basis = store.basis;
if ~isfield(store, "Hproj")
    problem_ = struct();
    problem_.M = M;
    problem_.hess = @(x, h, s) sylv_hess(x, h, A, B, C, s, struct("projehess", true));
    Hproj = hessianmatrix(problem_, X, basis); Hproj = symmetrize(Hproj);
    store.Hproj = Hproj;
end
if ~isfield(store, "Hcurv")
    problem_ = struct();
    problem_.M = M;
    problem_.hess = @(x, h, s) sylv_hess(x, h, A, B, C, s, struct("curvonly", true));
    Hcurv = hessianmatrix(problem_, X, basis); Hcurv = symmetrize(Hcurv);
    store.Hcurv = Hcurv;
end
if ~isfield(store, "H")
    dproj = eig(Hproj);
    [Vcurv,dcurv] = eig(Hcurv); dcurv = diag(dcurv);
    dhess = eig(Hproj + Hcurv);
    store.eigproj = dproj;
    store.eigcurv = dcurv;
    store.eighess = dhess;
    H = Hproj + Vcurv * diag(max(dcurv, zeros(size(dcurv)))) * Vcurv';
    store.H = H;
    
    % % Check if the eigenvalues of the curvature term are as in Proposition 2.30
    % [m, r] = size(X.U); n = size(X.V, 1); ell = size(store.AU, 3);
    % egrad.U = [reshape(bsxfun(@times, store.AU, diag(X.S)'), [m, r*ell]), -C.L];
    % egrad.V = [reshape(store.BV, [n, r*ell]), C.R]; perp_egrad = struct();
    % perp_egrad.U = egrad.U - X.U * [reshape(store.UtAUS, [r, r*ell]), -store.CLtU'];
    % perp_egrad.V = egrad.V - X.V * [reshape(store.VtBV, [r, r*ell]), store.CRtV'];
    % perp_egrad = perp_egrad.U * perp_egrad.V';
    % e = svd(perp_egrad); e = e(1:min([m-r, n-r]));
    % s = diag(X.S);
    % t = table2array(combinations(e,1./s));
    % eigen = t(:, 1) .* t(:, 2); eigen = sort(eigen, 'descend');
    % eigen = [-eigen; zeros(r*(m+n-r)-2*length(eigen), 1); flip(eigen)];
    % err = norm(eigen - dcurv) / norm(dcurv);
end
H = store.H;

% Get representation of eta
eta_vec = tangent2vec(M, X, basis, eta);

% Compute hessian approximation
Heta_vec = H * eta_vec;
Heta = lincomb(M, X, basis, Heta_vec);

end