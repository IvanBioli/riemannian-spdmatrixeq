function M_sym = symmetrize(M)
% symmetrize - Makes the symmetric matrix M "numerically symmetric" by
% computing M_sym = 0.5 * (M + M'). It also checks that M is "close to
% numerically symmetric".

    if length(size(M)) == 2 
        M_sym = 0.5 * (M + M');
    else
        M_sym = 0.5 * (pagetranspose(M) + M); 
    end
    norm_diff = norm(M_sym - M, 'fro');
    assert((norm_diff < 1e-10) || (norm_diff / norm(M, 'fro') < 1e-10))
end

