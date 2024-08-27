function [X_hat, frobnorm] = truncation(X, opts)
% truncation - Computes the truncation and the frobenius norm of (the
% matrix represented by) X
%
% Syntax:
%   [X_hat, frobnorm] = truncation(X, opts)
%
% Inputs:
%   - X   : Matrix to be truncated or structure representing the matrix to
%       be truncated, in factored form. The allowed factorizations are
%       X.U * X.V', X.U * X.S * X.V' and X.V * X.D * X.V'.
%   - opts: Options structure with field options.optionname, where
%       optionname is one of the following and the default value is indicated
%       between parentheses:
%           - method: Truncation method (see documentation of
%               truncation_rank_)
%           - tol: Truncation tolerance (see documentation of
%               truncation_rank_)
%           - rank_cap (Inf): Maximum rank of the truncated matrix
%           - SVD (false): Whether to return the truncation of X in SVD
%               format, if X is a struct representing X.U * X.V'
%
% Outputs:
%   - X_hat   : Structure representing the truncation of X
%   - frobnorm: Frobenius norm of X
%

method = opts.method;
tol = opts.tol;
if isfield(opts, 'rank_cap')
    rank_cap = opts.rank_cap;
else
    rank_cap = Inf;
end
frobnorm = NaN;
if isfield(X, 'U') && isfield(X, 'V') % Truncation of UV'
    [Q_U, R_U] = qr(full(X.U), 0);
    [Q_V, R_V] = qr(full(X.V), 0);
    [U_hat,S,V_hat] = svd(full(R_U * R_V'), 'econ'); S = diag(S);
    k_tilde = truncation_rank(S, method, tol);
    U_prime = Q_U * U_hat(:, 1:k_tilde);
    V_prime = Q_V * V_hat(:, 1:k_tilde);
    S_prime = diag(S(1:k_tilde));
    if isfield(opts, 'SVD') && opts.SVD
        X_hat.U = U_prime; X_hat.V = V_prime; X_hat.S = S_prime;
    else
        X_hat.U = U_prime * diag(sqrt(S(1:k_tilde)));
        X_hat.V = V_prime * diag(sqrt(S(1:k_tilde)));
    end
    frobnorm = sqrt(sum(S.^2));
elseif isfield(X, 'V') && isfield(X, 'D')% Truncation of VDV'
    [Q_V, R_V] = qr(full(X.V), 0);
    aux = full(R_V * X.D * R_V'); aux = 0.5 * (aux + aux');
    [U_hat,d_hat] = eig(aux); d_hat = diag(d_hat);
    frobnorm = sqrt(sum(d_hat.^2));
    % Sort eigenvalues by abs, i.e. in the same order as singular values
    [~, idx] = sort(abs(d_hat), 'descend');
    d_hat = d_hat(idx);
    U_hat = U_hat(:,idx);
    % Truncate
    k_tilde = truncation_rank(abs(d_hat), method, tol);
    X_hat.V = Q_V * U_hat(:, 1:k_tilde);
    X_hat.D = spdiags(d_hat(1:k_tilde), 0 , k_tilde, k_tilde);
elseif ~isstruct(X) % Truncation of a matrix
    [U_hat,S,V_hat] = svd(full(X), 'econ'); S = diag(S);
    k_tilde = truncation_rank(S, method, tol);
    X_hat=U_hat(:, 1:k_tilde) * diag(S(1:k_tilde)) * V_hat(:, 1:k_tilde)';
else
    error('Unsupported type of factorization to be truncated')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function k = truncation_rank(s, method, tol)
        % truncation_rank - Computes the truncation rank, based on the prescribed
        % truncation method and tolerance. Allows for multiple truncation methods
        % and tolerances, and then selects the minimum truncation rank among the
        % obtained ones. See documentation of truncation_rank_.
        
        k_vec = zeros(length(tol),1);
        for i = 1:length(tol)
            k_vec(i) = truncation_rank_(s, method(i), tol(i));
        end
        k = min(k_vec);
    end

    function k = truncation_rank_(s, method, tol)
        % truncation_rank_ - Computes the truncation rank based on the prescribed
        % method and tolerance. The truncation rank is always below rank_cap.
        %
        % Syntax:
        %   [k] = truncation_rank_(s, method, tol)
        %
        % Inputs:
        %   - s     : Array of singular values
        %   - method: Truncation method, i.e. criterion for selecting the
        %       truncation rank. It can be one of the following
        %           - "relF"/"absF": Relative/Absolute error in Frobenius norm
        %           - "rel2"/"abs2": Relative/Absolute error in 2-norm
        %           - "rank": Prescribed (maximum) truncation rank
        %   - tol   : Truncation tolerance
        %
        % Outputs:
        %   - k: Truncation rank
        %
        
        
        if method == "relF" % Relative error in Frobenius norm
            sq_err = cumsum([s(2:end).^2; 0], "reverse");
            sq_sum = sum(s.^2);
            sq_err = sq_err / sq_sum;
            k = find(sq_err <= tol^2, 1);
        elseif method == "absF" % Absolute error in Frobenius norm
            sq_err = cumsum([s(2:end).^2; 0], "reverse");
            k = find(sq_err <= tol^2, 1);
        elseif method == "rel2" % Relative error in 2-norm
            err = s / s(1);
            k = find(err <= tol, 1) - 1;
        elseif method == "abs2" % Absolute error in 2-norm
            err = s;
            k = find(err <= tol, 1) - 1;
        elseif method == "rank" % Maximum rank
            k = min([tol, size(X.U,2)]);
        else
            error("Unknown truncation method")
        end
        k = min([k, rank_cap, size(X.V,2)]);
    end

end

