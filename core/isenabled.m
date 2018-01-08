function rv = isenabled(mode, varargin)
    %   ISENABLED  Checks if mode exists in the cell-array varargin.
    %
    %   isenabled(mode,varargin{:}) return true or false.
    %   example:
    %
    %          varargin = {'Viz', 'ElementNumber', 'debug', [20,20]};
    %          isenabled('debug',varargin)
    %          ans =
    %               1
    %
    %   Author: Mirza Cenanovic (mirza.cenanovic@ju.se)
    %   Date: 2013-05-02
    if nargin < 1
        error('No arguments')
    end
    varargin = varargin{:};
    rv = 0;
    if any(strcmpi(varargin,mode))
        rv = 1;
    end
end