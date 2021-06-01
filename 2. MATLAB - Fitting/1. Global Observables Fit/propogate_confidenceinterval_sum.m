function [ds] = propogate_confidenceinterval_sum(dx)
%propogate_confidenceinterval_quotient computes confidence interval for a
%sum of all elements x.
%   Assumes elements in x are not correlated and follow a normal
%   distribution. Formula used then yields for the interval for the sum s =
%   sum x the following
%   ds = \sqrt( sum dx^2 )
%   Function allows for matrix input, where the summation are taken over
%   the second dimensions (i.e. rows).

ds = sqrt( sum( dx.^2, 2 ) );


% Other way to estimate error propogation (less precise)
%ds = sum(abs(dx),2);

end

