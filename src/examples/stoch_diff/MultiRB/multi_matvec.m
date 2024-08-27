function V = multi_matvec(K,M1,M2,v);

V=M2\v;
for k=2:length(K)
   V(:,k) = K{k}*V(:,1);
end
V=M1\V;
