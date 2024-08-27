function Tmn = TvecMat_sparse(m,n)
% TvecMat_sparse - Same as TvecMat, but returns a sparse matrix

d = m*n;

i = 1:d;
rI = 1+m.*(i-1)-(m*n-1).*floor((i-1)./n);
Tmn = sparse(rI, 1:d, ones(d, 1), d, d);
Tmn = Tmn';

