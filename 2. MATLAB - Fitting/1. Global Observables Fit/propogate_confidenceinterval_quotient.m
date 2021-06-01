function [dq] = propogate_confidenceinterval_quotient(q,x,dx,y,dy)
%propogate_confidenceinterval_quotient computes confidence interval for
%quotient q = x / y given confidence intervals (of same percentage) for x
%and y.
%   Assumes x and y are not correlated and follow a normal distribution,
%   and assumes q follows a normal distribution as well. Inputs are q:
%   mean value for the quotient (i.e. q = x / y), mean values x and y for
%   the variables and confidence interval lengths dx resp dy for x resp y.
%   Formula used is (std_q / q)^2 = (std_x / x)^2 + (std_y / y)^2. In case
%   x and y are correlated an additional term is present but that is not
%   programmed here.

dq = abs(q) .* sqrt( (dx ./ x).^2 + (dy ./ y).^2 );


% Other way to estimate error propogation (less precise)
%dq = abs(q) .* ( abs(dx ./ x) + abs(dy ./ y) );

end

