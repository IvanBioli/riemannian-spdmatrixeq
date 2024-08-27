function Y = stoch_matvec_m(X,G,K)
%STOCH_MATVEC efficient matvec for sum of kronecker product matrices
%   y = stoch_matvec(x,G,K);
%   input:
%          x    input vector of dimension (nKi $\times$ nGi)  
%          G    cell structure with nrv (nGi $\times$ nGi) matrices
%          K    cell structure with nrv (nKi $\times$ nKi) matrices 
%   output:
%          y    output vector of dimension (nKi $\times$ nGi)
%
% SIFISS function: DJS; 19 January 2013.
% Copyright (c) 2013 A. Bespalov, C.E. Powell, D.J. Silvester

dimk=length(K); dimg=length(G);
if dimk ~= dimg, error('incompatible cell dimensions'), end
% get dimensions
[n, dummy] = size(K{1});
[p, dummy] = size(G{1});
Y=zeros(size(X));
%
% loop over the number of matrices
for dim = 1:dimk
Y = Y + K{dim}*X*G{dim};
end
%
return
