function  y = Emergence_IO_BetaPDF(x,a,b,n)
% Way faster than the regular MATLAB "betapdf" function because:
%   - it does not check the inputs, but directly "repmat" the coefficients
%   - it does not have to check whether the inputs are smaller than 0
%   - it does not have to check the size of x

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