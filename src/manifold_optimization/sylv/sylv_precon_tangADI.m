function [xi, store] = sylv_precon_tangADI(X, eta, A, B, store, options, E, D)
% sylv_precon_tangADI - Applies steps of the tangADI algorithm
% to approximatively solve the projected generalized Sylvester equation
% Proj_X (A * xi * D + E * xi * B) = eta, for xi in Tang_X(M_r).
%
% Syntax:
%   [xi, store] = sylv_precon_tangADI(X, eta, A, B, store, options, E, D)
%
% Inputs:
%   - X      : Current iterate
%   - eta    : Vector of the tangent space at X to which the preconditioner
%              is applied
%   - A, B   : Coefficient matrices of the preconditioning operator
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure with fields p and q, which are vectors
%              containing the ADI shift parameters.
%   - D      : Coefficient matrix of the preconditioning operator. If not
%              specified it is assumed to be the identity matrix
%   - E      : Coefficient matrix of the preconditioning operator. If not
%              specified it is assumed to be the identity matrix
%
% Outputs:
%   - xi   : The result of the application of tangADI to eta
%   - store: Updated store struct
%

is_ED = nargin > 7;

% Definition of relevant quantities
U = X.U; V = X.V;
[m, r] = size(U); [n, ~] = size(V);
p = options.p; q = options.q; steps = length(p);

% Precomputations
AU = A * U; UtAU = U' * AU; UtAU = symmetrize(UtAU);
BV = B * V; VtBV = V' * BV; VtBV = symmetrize(VtBV);
if ~is_ED
    E = eye(m, 'like', A);
    D = eye(n, 'like', B);
    UtEU = eye(r, 'like', UtAU);
    VtDV = eye(r, 'like', VtBV);
else
    UtEU = U' * (E * U); UtEU = symmetrize(UtEU);
    VtDV = V' * (D * V); VtDV = symmetrize(VtDV);
end

% ADI iterations
for i = 1:steps
    % Compute the RHS
    if i == 1
        ZV = (p(i) - q(i)) * (eta.Up + X.U * eta.M);
        ZtU = (p(i) - q(i)) * (eta.Vp + X.V * eta.M');
        UtZV = (p(i) - q(i)) * eta.M;
    else
        Y = [X.U*xi.M + xi.Up, X.U]; W = [X.V, xi.Vp];
        ApY = (A - p(i)*E) * Y; BqW = (B + q(i)*D) * W;
        ZV = ApY * (BqW' * V);
        ZV = ZV - (q(i) - p(i)) * (eta.Up + X.U * eta.M);
        ZtU = BqW * (ApY' * U);
        ZtU = ZtU - (q(i) - p(i)) * (eta.Vp + X.V * eta.M');
        UtZV = U' * ZV;
    end
    
    % Update xi
    UtAU_p = UtAU - q(i)*UtEU; UtAU_p_d = decomposition(UtAU_p, 'chol');
    VtBV_q = VtBV + p(i)*VtDV; VtBV_q_d = decomposition(VtBV_q, 'chol');
    xi.Up = (A - q(i)*E) \ ZV / VtBV_q_d;
    xi.Up = xi.Up - U * (U' * xi.Up);
    xi.Vp = (B + p(i)*D) \ ZtU / UtAU_p_d;
    xi.Vp = xi.Vp - V * (V' * xi.Vp);
    xi.M = UtAU_p_d \ (UtZV - (U' * (A - q(i)*E) * xi.Up) * VtBV_q - ...
        ((V' * (B + p(i)*D) * xi.Vp) * UtAU_p)') / VtBV_q_d;
end

end