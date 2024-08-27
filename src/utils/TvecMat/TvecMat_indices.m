function rI = TvecMat_indices(m,n)
% TvecMat_indices - Same as TvecMat, but returns the indices corresponding
% to the permutation matrix in TvecMat.

d = m*n;

i = 1:d;
rI = 1+m.*(i-1)-(m*n-1).*floor((i-1)./n);
end
