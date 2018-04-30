function prod = norm_inner(a, b) 
% normalized inner product
%   a: ncell X t; b: t X n
%   Jinghao Lu 01/12/2016

    prod = a * b ./ (sqrt(sum(a .^ 2, 2)) * sqrt(sum(b .^ 2, 1)));
end
