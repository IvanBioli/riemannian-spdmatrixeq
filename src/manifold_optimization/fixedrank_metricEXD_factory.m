function M = fixedrank_metricEXD_factory(E, D, k, chol_E, chol_D)
% Manifold struct to optimize fixed-rank matrices w/ the geometry induced by P(X) = EXD
%
% function M = fixedrank_metricEXD_factory(m, n, k)
%
% Manifold of m-by-n real matrices of fixed rank k, embedded submanifold of
% the m-by-n real matrices with the metric induced by P(X) = EXD, where E
% and D are a symmetric positive definite matrices. See the Thesis for
% further details.
%
% A point X on the manifold is represented as a structure with three
% fields: U, S and V. The matrix U (mxk) has E-orthonormal columns, i.e.
% U'*E*U = I, the matrix and V (nxk) has D-orthonormal columns, while the
% matrix S (kxk) is any /diagonal/ full-rank matrix. Following the E,D-SVD
% formalism, X = U*S*V'. Note that the diagonal entries of S are not
% constrained to be nonnegative.
%
% Tangent vectors are represented as a structure with three fields: Up, M
% and Vp. The matrices Up (mxk) and Vp (mxk) obey Up'*E*U = 0 and V
% p'*D*V = 0. The matrix M (kxk) is arbitrary. Such a structure corresponds
% to the following tangent vector in the ambient space of mxn matrices:
%   Z = U*M*V' + Up*V' + U*Vp'
% where (U, S, V) is the current point and (Up, M, Vp) is the tangent
% vector at that point.
%
% Vectors in the ambient space are best represented as mxn matrices. If
% these are low-rank, they may also be represented as structures with
% U, S, V fields, such that Z = U*S*V'. There are no restrictions on what
% U, S and V are, as long as their product as indicated yields a real, mxn
% matrix.
%
% The chosen geometry yields a Riemannian submanifold of the embedding
% space R^(nxn) equipped with the inner product induced by P(X)=EXD (using
% the Frobenius metric as a base metric).
%

m = size(E, 1); n = size(D, 1);

M.name = @() sprintf(['Manifold of %dx%d matrices of rank %d, with ' ...
    'metric induced by P(X)=E*X*D for E,D symmetric positive ' ...
    'definite'], m, n, k);
M.E = E; M.D = D;
M.dim = @() (m+n-k)*k;

M.inner = @(x, d1, d2) inner(d1, d2);
    function inn = inner(d1, d2)
        inn = d1.M(:).'*d2.M(:);
        if isfield(d1, 'EUp')
            inn = inn + d1.EUp(:)' * d2.Up(:) + d1.DVp(:)' * d2.Vp(:);
        elseif isfield(d2, 'EUp')
            inn = inn + d1.Up(:)' * d2.EUp(:) + d1.Vp(:)' * d2.DVp(:);
        else
            d1_EUp = E * d1.Up; d1_DVp = D * d1.Vp;
            inn = inn + d1_EUp(:)' * d2.Up(:) + d1_DVp(:)' * d2.Vp(:);
        end
    end

M.norm = @(x, d) sqrt(inner(d, d));

M.dist = @(x, y) error('fixedrankembeddedfactory.dist not implemented yet.');

M.typicaldist = @() M.dim();

% Given Z in tangent vector format, projects the components Up and Vp
% such that they satisfy the tangent space constraints up to numerical
% errors. If Z was indeed a tangent vector at X, this should barely
% affect Z (it would not at all if we had infinite numerical accuracy).
M.tangent = @tangent;
    function Z = tangent(X, Z)
        Z.Up = Z.Up - X.U*(X.U'*(E * Z.Up));
        Z.Vp = Z.Vp - X.V*(X.V'*(D * Z.Vp));
        if isfield(Z, 'EUp')
            Z.EUp = E * Z.Up;
            Z.DVp = D * Z.Vp;
        end
    end

% For a given ambient vector Z, applies it to a matrix W. If Z is given
% as a matrix, this is straightforward. If Z is given as a structure
% with fields U, S, V such that Z = U*S*V', the product is executed
% efficiently.
    function ZW = apply_ambient(Z, W)
        if ~isstruct(Z)
            ZW = Z*W;
        else
            ZW = Z.U*(Z.S*(Z.V'*W));
        end
    end

% Same as apply_ambient, but applies Z' to W.
    function ZtW = apply_ambient_transpose(Z, W)
        if ~isstruct(Z)
            ZtW = Z'*W;
        else
            ZtW = Z.V*(Z.S'*(Z.U'*W));
        end
    end

% Orthogonal projection of an ambient vector Z represented as an mxn
% matrix or as a structure with fields U, S, V to the tangent space at
% X, in a tangent vector structure format.
M.proj = @projection;
    function Zproj = projection(X, Z)
        
        ZDV = apply_ambient(Z, X.DV);
        UtEZDV = X.EU'*ZDV;
        ZtEU = apply_ambient_transpose(Z, X.EU);
        
        Zproj.M = UtEZDV;
        Zproj.Up = ZDV  - X.U*UtEZDV;
        Zproj.Vp = ZtEU - X.V*UtEZDV';
        
    end

%     The two functions below need clarification about is "egrad" and "ehess",
%     i.e. what is considerd as "Euclidean". Formulas are provided in the
%     report, if need be to implement them.
%     M.egrad2rgrad = Not implemented yet.
%     M.ehess2rhess = Not implemented yet.

% Transforms a tangent vector Z represented as a structure (Up, M, Vp)
% into a structure with fields (U, S, V) that represents that same
% tangent vector in the ambient space of mxn matrices, as U*S*V'.
% This matrix is equal to X.U*Z.M*X.V' + Z.Up*X.V' + X.U*Z.Vp'. The
% latter is an mxn matrix, which could be too large to build
% explicitly, and this is why we return a low-rank representation
% instead. Note that there are no guarantees on U, S and V other than
% that USV' is the desired matrix. In particular, U and V are not (in
% general) orthonormal and S is not (in general) diagonal.
% (In this implementation, S is identity, but this might change.)
M.tangent2ambient_is_identity = false;
M.tangent2ambient = @tangent2ambient;
    function Zambient = tangent2ambient(X, Z)
        Zambient.U = [X.U*Z.M + Z.Up, X.U];
        Zambient.S = eye(2*k);
        Zambient.V = [X.V, Z.Vp];
    end

% This retraction is second order, following general results from
% Absil, Malick, "Projection-like retractions on matrix manifolds",
% SIAM J. Optim., 22 (2012), pp. 135-158.
%
% Notice that this retraction is only locally smooth, otherwise also
% the retraction of fixedrankembeddedfactory would be smooth. See
% fixedrankembeddedfactory.retraction for further details
M.retr = @retraction;
    function Y = retraction(X, Z, t)
        if nargin < 3
            t = 1.0;
        end
        
        % Mathematically, Z.Up is orthogonal to X.U, and likewise for
        % Z.Vp compared to X.V. Thus, in principle, we could call QR
        % on Z.Up and Z.Vp alone, which should be about 4 times faster
        % than the calls here where we orthonormalize twice as many
        % vectors. However, when Z.Up, Z.Vp are poorly conditioned,
        % orthonormalizing them can lead to loss of orthogonality
        % against X.U, X.V.
        [Qu, Ru] = qr_E([X.U, Z.Up]);
        [Qv, Rv] = qr_D([X.V, Z.Vp]);
        
        % Calling svds or svd should yield the same result, but BV
        % advocated svd is more robust, and it doesn't change the
        % asymptotic complexity to call svd then trim rather than call
        % svds. Also, apparently Matlab calls ARPACK in a suboptimal way
        % for svds in this scenario.
        % Notice that the parameter t appears only here. Thus, in princple,
        % we could make some savings for line-search procedures where we
        % retract the same vector multiple times, only with different
        % values of t. The asymptotic complexity remains the same though
        % (up to a constant factor) because of the matrix-matrix products
        % below which cost essentially the same as the QR factorizations.
        [U, S, V] = svd(Ru*[X.S + t*Z.M, t*eye(k); t*eye(k), zeros(k)]*Rv');
        
        Y.U = Qu*U(:, 1:k); Y.EU = E * Y.U;
        Y.V = Qv*V(:, 1:k); Y.DV = D * Y.V;
        Y.S = S(1:k, 1:k);
        
    end

% Less safe but much faster checksum, June 24, 2014.
% Older version right below.
M.hash = @(X) ['z' hashmd5([sum(X.U(:)) ; sum(X.S(:)); sum(X.V(:)) ])];
%M.hash = @(X) ['z' hashmd5([X.U(:) ; X.S(:) ; X.V(:)])];

M.rand = @random;
    function X = random()
        X.U = M.qr_E(randn(m, k)); X.EU = E * X.U;
        X.V = M.qr_D(randn(n, k)); X.DV = D * X.V;
        X.S = diag(sort(rand(k, 1), 1, 'descend'));
    end

% Generate a random tangent vector at X.
% Note: this may not be the uniform distribution over the set of
% unit-norm tangent vectors.
M.randvec = @randomvec;
    function Z = randomvec(X)
        Z.M  = randn(k);
        Z.Up = randn(m, k);
        Z.Vp = randn(n, k);
        Z = tangent(X, Z);
        nrm = M.norm(X, Z);
        Z.M  = Z.M  / nrm;
        Z.Up = Z.Up / nrm;
        Z.Vp = Z.Vp / nrm;
    end

M.lincomb = @lincomb;

M.zerovec = @(X) struct('M', zeros(k, k), 'Up', zeros(m, k), ...
    'Vp', zeros(n, k));

% New vector transport on June 24, 2014 (as indicated by Bart)
% Reference: Absil, Mahony, Sepulchre 2008 section 8.1.3:
% For Riemannian submanifolds of a Euclidean space, it is acceptable to
% transport simply by orthogonal projection of the tangent vector
% translated in the ambient space.
M.transp = @project_tangent;
    function Z2 = project_tangent(X1, X2, Z1)
        Z2 = projection(X2, tangent2ambient(X1, Z1));
    end

% It is sometimes useful to switch between representation of matrices
% as triplets or as full matrices of size m x n.
M.triplet2matrix = @triplet2matrix;
    function X_matrix = triplet2matrix(X_triplet)
        U = X_triplet.U;
        S = X_triplet.S;
        V = X_triplet.V;
        X_matrix = U*S*V';
    end

% Converts the SVD representation of a manifold's point to the ED-SVD
% representation, i.e. the one used for points on the manifold
M.SVD2repr = @SVD2repr;
    function Y = SVD2repr(X)
        [Qu, Ru] = qr_E([X.U]);
        [Qv, Rv] = qr_D([X.V]);
        [U, S, V] = svd(Ru * X.S * Rv');
        
        Y.U = Qu*U; Y.EU = E * Y.U;
        Y.V = Qv*V; Y.DV = D * Y.V;
        Y.S = S;
    end

% Computes the E-QR and the D-QR factorization, i.e. the QR
% factorization in the inner product induced by the matrices E, D.
% This can be done with the matrices E, D only, avoiding the Cholesky
% factorization (see details on Thesis). For simplicity, we still
% employ the Cholesky factorization.
if nargin < 4
    M.cholE = chol(E);
else
    M.cholE = chol_E;
end
if nargin < 5
    M.cholD = chol(D);
else
    M.cholD = chol_D;
end
M.qr_E = @qr_E;
M.qr_D = @qr_D;
    function [Q, R] = qr_E(Y)
        [Q, R] = qr(full(M.cholE * Y), 0);
        Q = M.cholE \ Q;
    end
    function [Q, R] = qr_D(Y)
        [Q, R] = qr(full(M.cholD * Y), 0);
        Q = M.cholD \ Q;
    end
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
if nargin == 3
    d.Up = a1*d1.Up;
    d.Vp = a1*d1.Vp;
    d.M  = a1*d1.M;
    if isfield(d1, 'EUp')
        d.EUp = a1*d1.EUp;
        d.DVp = a1*d1.DVp;
    end
elseif nargin == 5
    d.Up = a1*d1.Up + a2*d2.Up;
    d.Vp = a1*d1.Vp + a2*d2.Vp;
    d.M  = a1*d1.M  + a2*d2.M;
    if isfield(d1, 'EUp') && isfield(d2, 'EUp')
        d.EUp = a1*d1.EUp + a2*d2.EUp;
        d.DVp = a1*d1.DVp + a2*d2.DVp;
    end
else
    error('fixedrank.lincomb takes either 3 or 5 inputs.');
end

end
