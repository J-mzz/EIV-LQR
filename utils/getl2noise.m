function x = getl2noise(n, T, eps)
% GETL2NOISE  n x T noise whose columns satisfy ||x(:,k)||_2 <= eps.
%   A random Gaussian matrix is normalized by its spectral norm and scaled
%   by a random radius r <= eps, so every column has norm at most eps.

v = randn(n, T);
v_unit = v / norm(v, 2);
r = eps * (rand()^(1/n));
x = r * v_unit;

end