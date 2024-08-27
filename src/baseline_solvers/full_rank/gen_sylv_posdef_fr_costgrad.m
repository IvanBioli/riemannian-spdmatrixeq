function [f, g] = gen_sylv_posdef_fr_costgrad(X, A, B, C)
% gen_sylv_posdef_fr_costgrad - Computes the cost function for matrices in 
% full-rank format. The cost function is
%           f(X) = 0.5 * <calA(X), X> - <X, C>
% where calA(X) = A{1}*X*B{1} + ... + A{p}*X*B{p}

LX = sylv_op(A, B, X);
g = LX - C;
f = 0.5 * (X(:)' * LX(:)) - X(:)' * C(:);
end

