function [x, w, A] = gauss(n, a, b)

%------------------------------------------------------------------------------
% gauss.m
%------------------------------------------------------------------------------
%
% Purpose:
%
% Generates abscissas and weigths on I = [ a, b ] (for Gaussian quadrature).
%
%
% Syntax:
%
% [x, w, A] = gauss(n, a, b);
%
%
% Input:
%
% n    integer    Number of quadrature points.
% a    real       Left endpoint of interval.
% b    real       Right endpoint of interval.
%
%
% Output:
%
% x    real       Quadrature points.
% w    real       Weigths.
% A    real       Area of domain.
%------------------------------------------------------------------------------


% 3-term recurrence coefficients:
n = 1:(n - 1);
beta = 1 ./ sqrt(4 - 1 ./ (n .* n));

% Jacobi matrix:
J = diag(beta, 1) + diag(beta, -1); 


% Eigenvalue decomposition:

%
% e-values are used for abscissas, whereas e-vectors determine weights.
%

[V, D] = eig(J);
x = diag(D);


% Size of domain:
A = b - a;


% Change of interval:

%
% Initally we have I0 = [ -1, 1 ]; now one gets I = [ a, b ] instead.
%
% The abscissas are Legendre points.
%

if ~(a == -1 && b == 1)
  x = 0.5 * A * x + 0.5 * (b + a);
end


% Weigths:
w = V(1, :) .* V(1, :);
end