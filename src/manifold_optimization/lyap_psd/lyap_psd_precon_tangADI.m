function [xi, store] = lyap_psd_precon_tangADI(X, eta, A, store, options, M)
% lyap_psd_precon_tangADI - Applies steps of the tangADI algorithm
% to approximatively solve the projected generalized Sylvester equation
% Proj_X (A * xi * M + M * xi * A) = eta, for xi in Tang_X(Snr++).
%
% Syntax:
%   [xi, store] = lyap_psd_precon_tangADI(X, eta, A, store, options, M)
%
% Inputs:
%   - X      : Current iterate
%   - eta    : Vector of the tangent space at X to which the preconditioner
%              is applied
%   - A      : Coefficient matrix of the preconditioning operator
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure.  See mergeDefaultOptions
%   - M      : Coefficient matrix of the preconditioning operator. If not
%              specified it is assumed to be the identity matrix
%
% Outputs:
%   - xi   : The result of the application of tangADI to eta
%   - store: Updated store struct
%

is_M = nargin > 5;

% Definition of relevant quantities
V = X.V;
[n, r] = size(V);
p = options.p; steps = length(p);

% Precomputations
AV = A * V; VtAV = V' * AV; VtAV = symmetrize(VtAV);
if ~is_M
    M = eye(n, 'like', A);
    VtMV = eye(r, 'like', VtAV);
else
    VtMV = V' * (M * V); VtMV = symmetrize(VtMV);
end

% ADI iterations
for i = 1:steps
    % Compute the RHS
    if i == 1
        ZV = 2 * p(i) * (eta.Vp + V * eta.M);
        VtZV = 2 * p(i) * eta.M;
    else
        ApV = (A - p(i)*M) * V; ApVk = (A - p(i)*M) * xi.Vp;
        ZV = [ApV * xi.M + ApVk, ApV] * ([ApV, ApVk]' * V);
        ZV = ZV + 2 * p(i) * (eta.Vp + V * eta.M);
        VtZV = V' * ZV;
    end
    
    % Update xi
    VtAV_p = VtAV + p(i)*VtMV; VtAV_p_d = decomposition(VtAV_p, 'chol');
    xi.Vp = (A + p(i)*M) \ ZV / VtAV_p_d;
    xi.Vp = xi.Vp - V * (V' * xi.Vp);
    VtApVxi = V' * (A + p(i)*M) * xi.Vp;
    aux = VtApVxi * VtAV_p;
    xi.M = VtAV_p_d \ (VtZV - aux - aux') / VtAV_p_d;
end

end