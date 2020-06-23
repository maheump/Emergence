function pstar = Emergence_ProbWeighFun( p, gamma, p0 )
%EMERGENCE_PROBWEIGHFUN returns transformed probability values by a
%probability weighting function with parameters gamma and p0.
%   - "p": a scalar or an array specifying the raw probability value(s).
%   - "gamma": a scalar or an array specifying the slope of the function.
%   - "p0": a scalar or an array specifying the intercept of the function.
% See: https://en.wikipedia.org/wiki/Prospect_theory
% See: Gonzalez, R., & Wu, G. (1999). On the shape of the probability
%   weighting function. Cognitive psychology, 38(1), 129-166.
% 
% Copyright (c) 2020 Maxime Maheu

% Return transformed probability levels
pstar = logitinv(gamma .* logit(p) + (1 - gamma) .* logit(p0));

end

% Inverse of the logit transformation
function u = logitinv(v)
maxcut = -log(eps);
mincut = -log(1/realmin - 1);
u = 1 ./ (1 + exp(-max(min(v,maxcut),mincut)));
end

% Logit function
function a = logit(b)
a = log(b./(1-b));
a(real(a)~=a) = NaN;
end