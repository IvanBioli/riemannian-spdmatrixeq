function M = sympsdfixedrankembeddedfactory(n, k)
% Manifold struct to optimize fixed-rank matrices w/ an embedded geometry.
%
% function M = sympsdfixedrankembeddedfactory(n, k)
%
% Manifold of n-by-n real positive semidefinite matrices of fixed rank k.
% This follows the embedded geometry described in Bart Vandereycken's 2010
% PhD thesis: "Riemannian and multilevel optimization for rankconstrained
% matrix problems".
%
% A point X on the manifold is represented as a structure with two fields:
% V, D. The matrix V (nxk) has orthonormal colums, while the matrix D (kxk)
% is a diagonal matrix with positive entries. This represents X through its
% eigenvalue decomposition X = V*D*V'.
%
% Tangent vectors are represented as a structure with two fields: Vp and M.
% The matrix Vp (nxk) obeys Vp'*V = 0. The matrix M (kxk) is symmetric.
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
% space R^(nxn) equipped with the usual trace (Frobenius) inner product.
%

M.name = @() sprintf('Manifold of %dx%d symmetric positive semidefinite matrices of rank %d', n, n, k);

M.dim = @() n*k - k*(k-1)/2;

M.inner = @(x, d1, d2) d1.M(:).'*d2.M(:) + 2 * (d1.Vp(:).'*d2.Vp(:));

M.norm = @(x, d) sqrt(norm(d.M, 'fro')^2 + 2 * norm(d.Vp, 'fro')^2);

M.dist = @(x, y) error('sympsdfixedrankembeddedfactory.dist not implemented yet.');

M.typicaldist = @() M.dim();

% Given Z in tangent vector format, projects the components M and Vp
% such that they satisfy the tangent space constraints up to numerical
% errors. If Z was indeed a tangent vector at X, this should barely
% affect Z (it would not at all if we had infinite numerical accuracy).
M.tangent = @tangent;
    function Z = tangent(X, Z)
        Z.M = symm(Z.M);
        Z.Vp = Z.Vp - X.V*(X.V'*Z.Vp);
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
        symZV = apply_ambient(symZ, X.V);
        VtsymZV = X.V'*symZV;
        Zproj.M = VtsymZV;
        Zproj.Vp = symZV  - X.V*VtsymZV;
    end

M.egrad2rgrad = @projection;

% Given the Euclidean gradient at X and the Euclidean Hessian at X
% along H, where egrad and ehess are vectors in the ambient space and H
% is a tangent vector at X, returns the Riemannian Hessian at X along
% H, which is a tangent vector.
M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, H)
        
        % Euclidean part
        rhess = projection(X, ehess);
        
        % Curvature part
        T = apply_ambient(symm(egrad), H.Vp)/X.D;
        rhess.Vp = rhess.Vp + (T - X.V*(X.V'*T));
        
    end

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
        [Qv, Rv] = qr(full([X.V, Z.Vp]), 0);
        
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
        [~, idx] = sort(abs(d), 'descend');
        d = d(idx);
        V = V(:,idx);
        
        % Truncate to rank k
        Y.V = Qv*V(:, 1:k);
        Y.D = diag(d(1:k));
    end



% Less safe but much faster checksum, June 24, 2014.
% Older version right below.
M.hash = @(X) ['z' hashmd5([sum(X.V(:)); sum(X.D(:))])];
%M.hash = @(X) ['z' hashmd5([X.V(:) ; X.D(:))];

M.rand = @random;
% Factors V lives on Stiefel manifolds: reuse their random generator.
stiefeln = stiefelfactory(n, k);
    function X = random()
        X.V = stiefeln.rand();
        X.D = diag(rand(k, 1));
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


% It is sometimes useful to switch between representation of matrices
% as triplets or as full matrices of size nxn.
M.couple2matrix = @couple2matrix;
    function X_matrix = couple2matrix(X_couple)
        D = X_couple.D;
        V = X_couple.V;
        X_matrix = V*D*V';
    end

end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INUSL>

if nargin == 3
    d.Vp = a1*d1.Vp;
    d.M  = a1*d1.M;
elseif nargin == 5
    d.Vp = a1*d1.Vp + a2*d2.Vp;
    d.M  = a1*d1.M  + a2*d2.M;
else
    error('fixedrank.lincomb takes either 3 or 5 inputs.');
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
