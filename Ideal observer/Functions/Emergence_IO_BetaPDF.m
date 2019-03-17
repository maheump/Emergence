function  y = Emergence_IO_BetaPDF(x,a,b,n)
% EMERGENCE_BETAPDF returns the beta probability density as the MATLAB
% "betadf" function but in a much faster way because:
%   (1) it does not check the inputs, but directly "repmat" the coefficients
%   (2) it does not check whether the inputs are smaller than 0
%   (3) it does not check the size of x
%   - "x": values to be evaluated.
%   - "a": first parameter.
%   - "b": second parameter.
%   - "n": number of values to be evaluated (for efficiency).
% 
% Copyright (c) 2018 Maxime Maheu

if nargin < 4, n = numel(x); end

y = zeros(1, n);

a = repmat(a, [1, n]);
b = repmat(b, [1, n]);

y(a==1 & x==0) = b(a==1 & x==0);
y(b==1 & x==1) = a(b==1 & x==1);
y(a<1 & x==0) = Inf;
y(b<1 & x==1) = Inf;

k = a>0 & b>0 & x>0 & x<1;
a = a(k);
b = b(k);
x = x(k);

smallx = x<0.1;

loga = (a-1).*log(x);

logb = zeros(size(x), 'like', y);
logb( smallx) = (b( smallx)-1) .* log1p(-x( smallx));
logb(~smallx) = (b(~smallx)-1) .*  log(1-x(~smallx));

y(k) = exp(loga+logb - (gammaln(a)+gammaln(b)-gammaln(a+b)));

end