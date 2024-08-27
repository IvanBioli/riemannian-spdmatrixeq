function norm_est = hutch_normest(X, num_queries, hutch_args)
% stable_norm_fact - Estimates the frobenius norm of a factored matrix, 
% using the Hutch++ algorithm.
%
% Syntax:
%   [norm_est] = hutch_normest(X, type)
%
% Inputs:
%   - X             : Struct array representing the matrix in factored form. 
%       The allowed factorizations are X.U * X.V', X.U * X.S * X.V' and
%       X.V * X.D * X.V'.
%   - num_queries   : Total number of matrix-vector products to compute.
%   - hutch_args    : Cell array containing the options for Hutch++ (see
%       src/utils/randomizedLA/hutchplusplus.m).
%
% Outputs:
%   - norm_est      : Frobenius norm estimate
%

if nargin < 3
    hutch_args = {};
end

% Define the matVecOracle
[matVecOracleX, matVecOracleXt, dimensionX] = struct_to_matVecOracle(X);
   
% Call Hutch++ for estimating the trace of X'*X
matVecOracle =@(x) matVecOracleXt(matVecOracleX(x));
dimension = dimensionX(2);
norm_est = sqrt(hutchplusplus(matVecOracle, num_queries, dimension, hutch_args{:}));
end

