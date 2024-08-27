function inner = mat_inner(A, B)
% mat_inner - Computes the frobenius inner product of A and B.

inner = A(:)' * B(:);
end

