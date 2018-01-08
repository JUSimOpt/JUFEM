function [outA] = findValInVector(A, val, varargin)
% Find a given value in a vector.
% Reason: Removes the errors in search such as 48.000000001 by adding a tolerance. 

tol = 1e-14;
if nargin > 2
    tol = varargin{1};
end

outA = find(A >= val-tol & A <= val+tol); 
end

