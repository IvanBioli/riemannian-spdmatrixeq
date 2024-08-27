function [U, S, V, apxErr] = svdsketch_fun(A, tol, varargin)
%SVDSKETCH Matrix sketch in form of an SVD.
%
%   [U,S,V] = SVDSKETCH(A) returns the singular value decomposition (SVD)
%   of a low-rank matrix sketch of A. The matrix sketch only reflects the
%   most important features of A (up to a tolerance), which enables faster
%   calculation of the SVD of large matrices compared to using SVDS.
%
%   [U,S,V] = SVDSKETCH(A, tol) specifies a tolerance for the sketch of A
%   such that norm(U*S*V'-A,'fro')/norm(A,'fro') <= tol. The tol input must
%   satisfy sqrt(eps(class(A))) <= tol < 1, since SVDSKETCH does not detect
%   errors smaller than sqrt(eps(class(A))). The default is
%   eps(class(A))^(1/4).
%
%   [U,S,V] = SVDSKETCH(A, tol, Name, Value) configures additional optional  
%   options specified by one or more name-value pair arguments:
%
%   'MaxSubspaceDimension' - Maximum subspace dimension (maxDim) for the
%                            algorithm. Specified as a positive integer to
%                            control the memory consumption of the
%                            algorithm. The default value is min(size(A)).
%              'BlockSize' - Block size of the algorithm. A larger block
%                            size can minimize the number of needed
%                            iterations, but might also add more
%                            information than necessary to achieve
%                            convergence. Must be a positive integer. The
%                            default value is
%                            min(max(floor(0.1*m), 5), maxDim), where
%                            m = size(A,1).
%          'MaxIterations' - Maximum number of iterations for the
%                            algorithm. More iterations can produce a
%                            higher quality matrix sketch at the cost of
%                            more execution time. Must be a positive
%                            integer. The default value is
%                            MaxSubspaceDimension divided by BlockSize.
%     'NumPowerIterations' - Number of power iterations performed within
%                            each iteration of the algorithm. Power
%                            iterations improve the orthogonality of the U
%                            and V outputs. Must be a nonnegative integer,
%                            but is typically 0, 1, or 2, as larger values
%                            can adversely contribute to round-off error.
%                            The default value is 1.
%
%   [U,S,V,apxErr] = SVDSKETCH(...) additionally returns the vector apxErr,
%   whose entries represent the relative approximation error in each
%   iteration, norm(U*S*V'-A,'fro')/norm(A,'fro'). The length of apxErr is
%   equal to the number of iterations, and apxErr(end) is the relative
%   approximation error of the output.
%
%   See also: svds, eigs, svd, pca
%
%   Note: PCA is in the Statistic and Machine Learning Toolbox

% Copyright 2019-2020 The MathWorks, Inc.

% REFERENCES:
%   W. Yu, Y. Gu, and Y. Li, Efficient randomized algorithms for
%   the fixed-precision low-rank approximation, SIAM J. Matrix Anal.
%   Appl., 39(3), 2018, pp. 1339–1359.

% Check provided tolerance or set to default if not provided
if nargin>1
    tol = full(tol);
    if ~isnumeric(tol) || ~isreal(tol) ||~isscalar(tol) || ...
            tol < sqrt(eps(class(A))) || ~(tol < 1)
        error(message('MATLAB:svdsketch:InvalidTolerance'));
    end
else
    tol = eps(class(A))^(1/4);
end

% Parse optional inputs and set default values if values are not provided
[bSize, maxIter, maxDim, p] = parseInputs(size(A), varargin);

% The actual work is done by computing a Q*B approximation of A such that B
% is much smaller than A.
[Q, B, apxErr] = randsvdImpl(A, tol, bSize, maxIter, maxDim, p);

% B can be sparse for some sparse inputs, e.g., A = sparse([0,1;1,0])
% Compute the SVD of the smaller matrix B, where Q*B approximates A and
% satisfies the approximation error tolerance "tol"
[U, S, V] = svd(full(B), 'econ');
U = Q*U;
end

function [Q, B, apxErr] = randsvdImpl(A, tol, bSize, maxIter, maxDim, p)

classA = class(A);
% Set randstream for reproducibility
randStr = RandStream('threefry', 'Seed', 0);

bSize0 = bSize;
[m, n] = size(A);
Q = zeros(m, 0, classA);
B = zeros(0, n, classA);
apxErr = zeros(0, 1, classA);

normA = norm(A, 'fro');
normB = 0;
% If normA is zero, no computation needs to be performed and we return
% early
if normA==0
    return
end

% Compute QB-factorization of A
for i = 1:maxIter
    Omega_i = randn(randStr, n, bSize, classA);
    Qi = orth(A*Omega_i - Q*(B*Omega_i));
    
    % Power iterations to increase accuracy
    for j = 1:p
        Qi = orth(A'*Qi - B'*(Q'*Qi));
        Qi = orth(A*Qi - Q*(B*Qi));
    end
    
    % Orthogonalize Qi and compute Bi
    Qi = orth(Qi - Q*(Q'*Qi));
    Bi = Qi'*A;
    
    % Append blocks to Q and B
    Q = [Q, Qi]; %#ok<AGROW>
    B = [B; Bi]; %#ok<AGROW>
    bDim = size(B, 1);
    
    % Accumulate relative approximation error in the Frobenius norm:
    % apxErr = ||A - Q*B||/||A|| using the following equalities that are
    % derived in above mentioned reference:
    % ||A - Q*B||^2 = ||A||^2 - ||B||^2             (*)
    %               = ||A||^2 - sum(||Bi||^2)       (**)
    % First, we sum up the ||Bi||^2 in (**) via HYPOT to minimize overflow
    normB = hypot(normB, norm(Bi, 'fro'));
    % Secondly, we use
    % ||A||^2 - ||B||^2 = (||A|| - ||B||)*(||A|| + ||B||),
    % of which we calculate the square root and normalize by ||A||.
    apxErr(i,1) = sqrt(abs(normA-normB)*(normA+normB))/normA;
    
    % Check for convergence of relative approximation error
    if apxErr(i) <= tol
        break
    end
    
    % Sometimes the last step moves us away from the best achievable
    % result, don't go that step.
    % This should not be needed often, as it's intended to handle edge
    % cases, where we experience unfortunate round-off error accumulations.
    if (i > 1) && (apxErr(i) > apxErr(i-1))
        Q(:, end-bSize+1:end) = [];
        B(end-bSize+1:end, :) = [];
        apxErr(i) = [];
        break
    end
    
    % Increase the block size, incrementing by the original block size,
    % bSize0, to speed up convergence if apxErr does not decrease by half.
    if (i > 1) && (apxErr(i) > apxErr(i-1)/2)
        bSize = min(bSize + bSize0, max(maxDim - bDim, 1));
    end
    
    % Make sure maxDim is not exceeded on last iteration
    if bDim + bSize > maxDim
        bSize = maxDim - bDim;
    end
    
    % If bSize is set to 0, break
    if bSize == 0
        break
    end
end
end

function [bSize, maxIter, maxDim, p] = parseInputs(sizeA, args)

% Set default values for maxDim and p
maxDim = max(min(sizeA), 1);
p = 1;
% Set maxIter and bSize to [] for now to ensure user provided value does
% not get overwritten by updated derived value later.
maxIter = [];
bSize = [];

% Parse NVPs if provided
if ~isempty(args)
    nameList = {'BlockSize', 'MaxIterations', 'MaxSubspaceDimension', ...
        'NumPowerIterations'};
    
    if mod(length(args), 2) ~= 0
        error(message('MATLAB:svdsketch:NameWithoutValue'));
    end
    
    for i=1:2:(length(args))
        indName = strncmpiWithInputCheck(args{i}, nameList);
        
        if nnz(indName) ~= 1
            % Let validatestring give the error message
            validatestring(args{i}, nameList);
        else
            indName = find(indName);
        end
        
        value = args{i+1};
        
        switch indName
            case 1 % 'BlockSize'
                if ~isScalarInt(value) || ~(value > 0)
                    error(message('MATLAB:svdsketch:InvalidBsize'));
                end
                bSize = double(full(value));
            case 2 % 'MaxIterations'
                if ~isScalarInt(value) || ~(value > 0)
                    error(message('MATLAB:svdsketch:InvalidMaxiter'));
                end
                maxIter = double(full(value));
            case 3 % 'MaxSubspaceDimension'
                if ~isScalarInt(value) || ~(value > 0)
                    error(message('MATLAB:svdsketch:InvalidMaxdim'));
                end
                maxDim = double(full(value));
            case 4 % 'NumPowerIterations'
                if ~isScalarInt(value) || ~(value >= 0)
                    error(message('MATLAB:svdsketch:InvalidP'));
                end
                p = double(full(value));
        end
    end
end

% Compute bSize if not provided as NVP
if isempty(bSize)
    bSize = min(max(floor(0.1*sizeA(1)), 5), maxDim);
end
% Compute maxIter if not provided as NVP
if isempty(maxIter)
    maxIter = ceil(maxDim/bSize);
end
end

function index = strncmpiWithInputCheck(arg, list)
% Check that the input is either a row char vector or a scalar string. Use
% that valid input to compare it against a provided list and return the
%index. Return [] if validation failed.
index = [];
if (ischar(arg) && isrow(arg)) || (isstring(arg) && isscalar(arg))
    index = startsWith(list, arg, 'IgnoreCase', true);
end
end

function tf = isScalarInt(x)
% Checks for scalar integers.
tf = isscalar(x) && isnumeric(x) && isreal(x) && fix(x) == x;
end
