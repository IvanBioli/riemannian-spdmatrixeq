function Y = adi_lyap(A, X, p, options, M, N)
% adi_lyap - Applies the ADI preconditioner, possibly bilinear to X, i.e. 
% applies ADI to approximately solve.
%       A*X*M' + M*X*A' - N{1}*X*N{1}' - N{p}*X*N{p}'  = X.V*X.D*X.V'
% If N is present, the bilinear ADI iteration from Benner and Breiten,
% 2013, is applied. Otherwise, standard fADI is employed. Truncation can be
% employed in the case of bilinear ADI iteration to avoid excessive
% intermediate ranks growth.
%
% Syntax:
%   [Y] = adi_lyap(A, X, p, options, M, N)
%   [Y] = adi_lyap(A, X, p, options, M)
%   [Y] = adi_lyap(A, X, p, options)
%   [Y] = adi_lyap(A, X, p)
%
% Inputs:
%   - A      : Coefficient matrix
%   - X      : Struct of the RHS
%   - p      : ADI shifts
%   - options: Options structure
%   - M      : Coefficient matrix
%   - N      : Cell array containing coefficient matrices of the bilinear
%              terms
%
% Outputs:
%   - Y: Approximate solution in factored form  V * D * V'
%

if nargin < 4
    options = struct();
end
options = mergeDefaultOptions(options);

% Get dimensions and parse arguments
[n, r0] = size(X.V);
is_M = (nargin > 4) && ~isempty(M);
if ~is_M
    M = speye(n);
end

if nargin < 6 % Standard fADI iteration
    if is_M
        Y.V = fadi_lyap(A, X.V, p, M);
    else
        Y.V = fadi_lyap(A, X.V, p);
    end
    Y.D = kron(speye(length(p)), X.D);
else
    n_terms = length(N);
    n_iters = length(p); % Number of ADI iterations
    
    % Initialize Y = 0 and perform first iteration
    Y.V = sqrt(2 * p(1)) * ((A + p(1) * M) \ X.V);
    Y.D = X.D;
    
    for i = 2:n_iters
        % Preallocating space for Vnext
        r = size(Y.V, 2);
        Vnext = zeros(n, (n_terms + 1) * r + r0);
        
        % ADI iteration
        % Compute Vnext
        Vnext(:, 1:r) = (A - p(i) * M) * Y.V;
        for j = 1:n_terms
            Vnext(:, j * r + 1 : (j + 1) * r) =  sqrt(2 * p(i)) * (N{j} * Y.V); 
        end
        Vnext(:, (n_terms + 1) * r + 1:end) = sqrt(2 * p(i)) * X.V;
        Vnext = (A + p(i) * M) \ Vnext;
        Y.V = Vnext;
        % Compute Dnext
        dim_Dnext = (n_terms + 1) * size(Y.D, 2) + r0;
        d_next = [diag(Y.D); repmat(diag(Y.D), n_terms, 1); diag(X.D)];
        Dnext = spdiags(d_next, 0, dim_Dnext, dim_Dnext);
        Y.D = Dnext;
    
        % Possibly truncate
        if options.trunc_tol
            truncation_ops.rel = true;
            truncation_ops.tol = options.trunc_tol;
            Y = truncation(Y, truncation_ops);
        end
    end
end
end

function options = mergeDefaultOptions(options)
% mergeDefaultOptions - Merges the option structure given by the user with 
%   the default option structure. User-defined options have precedence 
%   over the default ones. 
%
% Inputs:
%   - options: Options structure with fields options.optionname, where
%       optionname is one of the following and the default value is indicated
%       between parentheses:
%           - trunc_tol (false): Relative truncation tolerance. False if no
%             truncation is applied.
%
% Outputs:
%   - options: Options structure containing all the fields above, with
%     default values in fields that have not been provided by the user
%


    default_options = struct();
    default_options.trunc_tol = false;
    % Merge options with defaults
    options = mergeOptions(default_options, options);
end