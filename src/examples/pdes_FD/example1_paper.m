function [data, u] = example1_paper(a, n_terms)
% example1 - Returns the data (for assemble_pdes_FD.m) and the boundary 
% data for example 1 in numerical experiments. The diffusion, source and
% boundary data are:
%       k(x,y) = 1 + sum_{i=1}^{n_terms} a^i / i! * x^i * y^i
%       f(x,y) = 0
%       g(x,y) = ??

data.eps = 1;
data.k_coeff = a.^(1:n_terms) ./ factorial(1:n_terms);
for i = 1:n_terms
    data.k_x{i} =@(x) x.^i;
    data.k_y{i} =@(y) y.^i;
end
data.k =@(x,y) truncated_exp(x, y, a, n_terms);
data.k_mean = data.eps + sum(data.k_coeff ./ ((2:n_terms+1).^2), 'all');

u =@(x,y) nan;
% data.g =@(x,y) cos(x) .* sin(y);
data.g =@(x,y) exp(-a * (x+1) .* y);
% data.g =@(x,y) exp(-a * x.* y);

% Dominant part
data_dominant = struct();
if a > 1
    idx = length(data.k_coeff);
else
    idx = 1;
end
data_dominant.eps = 0;
data_dominant.k_coeff = [1];
data_dominant.k_x = cell(1); data_dominant.k_x{1} =@(x) data.eps + sqrt(data.k_coeff(idx)) * data.k_x{idx}(x);
data_dominant.k_y = cell(1); data_dominant.k_y{1} =@(x) data.eps + sqrt(data.k_coeff(idx)) *data.k_y{idx}(x);
data_dominant.k =@(x, y) data_dominant.k_x{1}(x) .* data_dominant.k_y{1}(y);
data.dominant = data_dominant;
end

function r = truncated_exp(x, y, a, k)
% truncated_exp - Evaluation of the degree-k Taylor expansion of exp(a*x*y) 

r = 1;
for i = 1:k
    r = r + a^i/factorial(i) * x.^i .* y.^i;
end
end

