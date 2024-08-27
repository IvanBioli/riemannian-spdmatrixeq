function M = sympsdfixedrank_metricMXM_factory (Mass, k)
% Manifold struct to optimize fixed-rank matrices w/ the geometry induced by P(X) = Mass*X*Mass.
%
% function M = sympsdfixedrankembeddedfactory(n, k)
%
% Manifold of n-by-n real positive semidefinite matrices of fixed rank k,
% as an embedded submanifold of the n-by-n real matrices with the metric
% induced by P(X) = Mass*X*Mass, where Mass is a symmetric positive
% definite matrix. See the Thesis for further details.
%
% A point X on the manifold is represented as a structure with two fields:
% V, D. The matrix V (nxk) has Mass-orthonormal colums, i.e. V'*Mass*V = I,
% while the matrix D (kxk) is a diagonal matrix with positive entries. This
% represents X through the decomposition X = V*D*V'. Note that this is not
% the eigenvalue decomposition of X.
%
% Tangent vectors are represented as a structure with two fields: Vp and M.
% The matrix Vp (nxk) obeys Vp'*Mass*V = 0. The matrix M (kxk) is symmetric.
% Such a structure corresponds to the following tangent vector in the
% ambient space of nxn matrices:
%           Z = V*M*V' + Vp*V' + V*Vp'
% where (V, D) is the current point and (M, Vp) is the tangent
% vector at that point.
%
% Vectors in the ambient space are best represented as nxn matrices. If
% these are low-rank, they may also be represented as structures with
% U, S, V fields, such that Z = U*S*V'. There are no restrictions on what
% U, S and V are, as long as their product as indicated yields a real, nxn
% matrix.
%
% The chosen geometry yields a Riemannian submanifold of the embedding
% space R^(nxn) equipped with the inner product induced by P(X)=Mass*X*Mass
% (using the Frobenius metric as a base metric).
%

n = size(Mass, 1);
M.Mass = Mass;
M.name = @() sprintf(['Manifold of %dx%d symmetric positive ...' ...
    'semidefinite matrices of rank %d with metric induced by ...' ...
    'P(X)=Mass*X*Mass for Mass symmetric positive'], n, n, k);
M.dim = @() n*k - k*(k-1)/2;

M.inner = @(x, d1, d2) inner(d1, d2);
    function inn = inner(d1, d2)
        inn = d1.M(:).'*d2.M(:);
        if isfield(d1, 'MVp')
            inn = inn + 2 * d1.MVp(:)' * d2.Vp(:);
        elseif isfield(d2, 'MVp')
            inn = inn + 2 * d1.Vp(:)' * d2.MVp(:);
        else
            d1_MVp = Mass * d1.Vp;
            inn = inn + 2 * d1_MVp(:)' * d2.Vp(:);
        end
    end

M.norm = @(x, d) sqrt(inner(d, d));

M.dist = @(x, y) error('sympsdfixedrank_metricMXM_factory.dist not implemented yet.');

M.typicaldist = @() M.dim();

% Given Z in tangent vector format, projects the components M and Vp
% such that they satisfy the tangent space constraints up to numerical
% errors. If Z was indeed a tangent vector at X, this should barely
% affect Z (it would not at all if we had infinite numerical accuracy).
M.tangent = @tangent;
    function Z = tangent(X, Z)
        Z.M = symm(Z.M);
        Z.Vp = Z.Vp - X.V*(X.V'*(Mass * Z.Vp));
        if isfield(Z, 'MVp')
            Z.MVp = Mass * Z.Vp;
        end
    end

% For a given ambient vector Z, applies it to a matrix W. If Z is given
% as a matrix, this is straightforward. If Z is given as a structure
% with fields U, S, V such that Z = U*S*V', or as a structure with
% fields V, D such that Z = V*D*V' with symmetric D, the product is
% executed efficiently.
    function ZW = apply_ambient(Z, W)
        if ~isstruct(Z)
            ZW = Z*W;
        else
            if isfield(Z, 'S')
                ZW = Z.U*(Z.S*(Z.V'*W));
            elseif isfield(Z, 'D')
                ZW = Z.V*(Z.D*(Z.V'*W));
            else
                error('Unrecognized struct for representing ambient space matrices.')
            end
        end
    end

% Same as apply_ambient, but applies Z' to W.
    function ZtW = apply_ambient_transpose(Z, W)
        if ~isstruct(Z)
            ZtW = Z'*W;
        else
            if isfield(Z, 'S')
                ZtW = Z.V*(Z.S'*(Z.U'*W));
            elseif isfield(Z, 'D')
                ZtW = Z.V*(Z.D*(Z.V'*W));
            else
                error('Unrecognized struct for representing ambient space matrices.')
            end
        end
    end

% Orthogonal projection of an ambient vector Z represented as an nxn
% matrix or as a structure with fields U, S, V (or V, D) to the tangent
% space at X, in a tangent vector structure format.
M.proj = @projection;
    function Zproj = projection(X, Z)
        symZ = symm(Z);
        symZMV = apply_ambient(symZ, X.MV);
        VtMsymZMV = X.MV'*symZMV;
        Zproj.M = VtMsymZMV;
        Zproj.Vp = symZMV  - X.V*VtMsymZMV;
    end

%     The two functions below need clarification about is "egrad" and "ehess",
%     i.e. what is considerd as "Euclidean". Formulas are provided in the
%     report, if need be to implement them.
%     M.egrad2rgrad = Not implemented yet.
%     M.ehess2rhess = Not implemented yet.


% Transforms a tangent vector Z represented as a structure (M, Vp)
% into a structure with fields (D, V) that represents that same
% tangent vector in the ambient space of nxn matrices, as V*D*V'.
% This matrix is equal to X.V*Z.M*X.V' + Z.Vp*X.V' + X.V*Z.Vp'. The
% latter is an nxn matrix, which could be too large to build
% explicitly, and this is why we return a low-rank representation
% instead. Note that there are no guarantees on D and V other than
% that VDV' is the desired matrix and D=D'. In particular, V is not (in
% general) orthonormal and D is not (in general) diagonal.
M.tangent2ambient_is_identity = false;
M.tangent2ambient = @tangent2ambient;
    function Zambient = tangent2ambient(X, Z)
        Zambient.V = [X.V, Z.Vp];
        I = eye(k); O = zeros(k);
        Zambient.D = [Z.M, I;...
            I,   O];
    end

% This retraction is second order, following general results from
% Absil, Malick, "Projection-like retractions on matrix manifolds",
% SIAM J. Optim., 22 (2012), pp. 135-158.
M.retr = @retraction;
    function Y = retraction(X, Z, t)
        if nargin < 3
            t = 1.0;
        end
        
        % Mathematically, Z.Vp is orthogonal to X.V. Thus, in principle, we
        % could call QR on Z.Up and Z.Vp alone, which should be about 2
        % times faster than the calls here where we orthonormalize twice as
        % many vectors. However, when Z.Vp is poorly conditioned,
        % orthonormalizing them it lead to loss of orthogonality against
        % X.V.
        [Qv, Rv] = qr_M([X.V, Z.Vp]);
        
        % Matematically, C is already symmetric. However, to make sure that
        % this is used when computing its eigenvalue decomposition, we make
        % sure it is also numerically symmetric
        C = Rv*[X.D + t*Z.M, t*eye(k); t*eye(k), zeros(k)]*Rv';
        C = symm(C);
        
        % Notice that the parameter t appears only here. Thus, in princple,
        % we could make some savings for line-search procedures where we
        % retract the same vector multiple times, only with different
        % values of t. The asymptotic complexity remains the same though
        % (up to a constant factor) because of the matrix-matrix products
        % below which cost essentially the same as the QR factorizations.
        [V,d] = eig(C); d = diag(d);
        % Sort eigenvalues by abs, i.e. same order as singular values
        [~, idx] = sort(d, 'descend');
        d = d(idx);
        V = V(:,idx);
        
        % Truncate to rank k
        assert(d(k) > 0)
        Y.V = Qv*V(:, 1:k); Y.MV = Mass * Y.V;
        Y.D = diag(d(1:k));
    end



% Less safe but much faster checksum, June 24, 2014.
% Older version right below.
M.hash = @(X) ['z' hashmd5([sum(X.V(:)); sum(X.D(:))])];
%M.hash = @(X) ['z' hashmd5([X.V(:) ; X.D(:))];

M.rand = @random;
    function X = random()
        X.V = M.qr_M(randn(n, k)); X.MV = Mass * X.V;
        X.D = diag(sort(rand(k, 1), 1, 'descend'));
    end

% Generate a random tangent vector at X.
% Note: this may not be the uniform distribution over the set of
% unit-norm tangent vectors.
M.randvec = @randomvec;
    function Z = randomvec(X)
        Z.M  = randn(k);
        Z.Vp = randn(n, k);
        Z = tangent(X, Z);
        nrm = M.norm(X, Z);
        Z.M  = Z.M  / nrm;
        Z.Vp = Z.Vp / nrm;
    end

M.lincomb = @lincomb;

M.zerovec = @(X) struct('M', zeros(k, k), 'Vp', zeros(n, k));

% New vector transport on June 24, 2014 (as indicated by Bart)
% Reference: Absil, Mahony, Sepulchre 2008 section 8.1.3:
% For Riemannian submanifolds of a Euclidean space, it is acceptable to
% transport simply by orthogonal projection of the tangent vector
% translated in the ambient space.
M.transp = @project_tangent;
    function Z2 = project_tangent(X1, X2, Z1)
        Z2 = projection(X2, tangent2ambient(X1, Z1));
    end

% Converts the eigen representation of a manifold's point to the
% M-orthgonal eigen representation, i.e. the one used for points on
% the manifold
M.eig2repr = @eig2repr;
    function Y = eig2repr(X)
        [Qv, Rv] = qr_M([X.V]);
        [V, D] = eig(symm(Rv * X.D * Rv'));
        
        Y.V = Qv*V; Y.MV = Mass * Y.V;
        Y.D = D;
    end


% It is sometimes useful to switch between representation of matrices
% as couples or as full matrices of size n x n.
M.couple2matrix = @couple2matrix;
    function X_matrix = couple2matrix(X_couple)
        D = X_couple.D;
        V = X_couple.V;
        X_matrix = V*D*V';
    end

% Computes the Mass-QR factorization, i.e. the QR factorization in the
% inner product induced by the Mass matrix.
% This can be done with the matrix Mass only, avoiding its Cholesky
% factorization (see details on Thesis). For simplicity, we still
% employ the Cholesky factorization
M.cholMass = chol(Mass);
M.qr_M = @qr_M;
    function [Q, R] = qr_M(Y)
        [Q, R] = qr(full(M.cholMass * Y), 0);
        Q = M.cholMass \ Q;
    end
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>
if nargin == 3
    d.Vp = a1*d1.Vp;
    d.M  = a1*d1.M;
    if isfield(d1, 'MVp')
        d.MVp = a1*d1.MVp;
    end
elseif nargin == 5
    d.Vp = a1*d1.Vp + a2*d2.Vp;
    d.M  = a1*d1.M  + a2*d2.M;
    if isfield(d1, 'MVp') && isfield(d2, 'MVp')
        d.MVp = a1*d1.MVp + a2*d2.MVp;
    end
else
    error('sympsdfixedrank_metricMXM_factory.lincomb takes either 3 or 5 inputs.');
end

end

% Symmetrization operator
function symX = symm(X)
if ~isstruct(X)
    symX = .5*(X+X');
else
    if isfield(X, 'S')
        symX.U = [X.U, X.V];
        S = .5*X.S; symX.S = blkdiag(S, S');
        symX.V = [X.V, X.U];
    elseif isfield(X, 'D')
        symX = X;
    else
        error('Unrecognized struct for representing ambient space matrices.')
    end
end
end
