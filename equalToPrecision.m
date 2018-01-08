function ind = equalToPrecision(v,val,varargin)

tol = 1e-12;
if nargin > 2
    tol = varargin{1};
end

ind =  v < val+tol & v > val-tol;

end