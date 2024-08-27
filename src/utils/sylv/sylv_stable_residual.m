function res = sylv_stable_residual(A, B, C, X)
% sylv_stable_residual - Computes the Frobenius norm of the residual
%       R = A{1} * X * B{1} + ... + A{p} * X * B{p} - C.L * C.R',
% in a numerically stable way.
%
% Syntax:
%   [res] = sylv_stable_residual(A, B, C, X)
%
% Inputs:
%   - A, B: Coefficient matrices
%   - C   : Structure representing the RHS F = C.L * C.R'
%   - X   : Structure representing the current iterate in factored form.
%       Allowed factorizations are X.U * X.S * X.V' and X.U * X.V'.
%
% Outputs:
%   - res: Residual norm
%

X_UV = struct();
if isfield(X, 'S')
    X_UV.U = X.U * X.S;
    X_UV.V = X.V;
else
    X_UV = X;
end
LX = sylv_op_lr(A, B, X_UV);
LX_C.U = [LX.U, C.L]; LX_C.V = [LX.V, -C.R];
res = stable_norm_fact(LX_C);
end

