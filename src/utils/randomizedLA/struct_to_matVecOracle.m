function [matVecOracleX, matVecOracleXt, dimensionX] = struct_to_matVecOracle(X)
% Defining the mat-vec oracle and the dimension

if isfield(X, 'V') && ~isfield(X, 'U')
    matVecOracleX =@(x) X.V * (X.D * (X.V' * x));
    matVecOracleXt =@(x) X.V * (X.D' * (X.V' * x));
    dimensionX = [size(X.V, 1) size(X.V, 1)];
elseif isfield(X, 'U') && isfield(X, 'S')
    matVecOracleX =@(x) X.U * (X.S * (X.V' * x));
    matVecOracleXt =@(x) X.V * (X.S' * (X.U' * x));
    dimensionX = [size(X.U, 1) size(X.V, 1)];
elseif isfield(X, 'U') && ~isfield(X, 'S')
    matVecOracleX =@(x) X.U * (X.V' * x);
    matVecOracleXt =@(x) X.V * (X.U' * x);
    dimensionX = [size(X.U, 1) size(X.V, 1)];
elseif ~isstruct(X)
    matVecOracleX =@(x) X * x;
    matVecOracleXt =@(x) X' * x;
    dimensionX = size(X);
else
    error("Unrecognized struct for X")
end

end

