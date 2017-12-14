function r = roundNR(in, dec)
% round up numbers
% in, as a scalar, vector or matrix
% dec, number of decimals

r = round((in.*10^dec)/10^dec,dec);


end
