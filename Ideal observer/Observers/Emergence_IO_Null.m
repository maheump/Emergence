function [ pY, pAgY ] = Emergence_IO_Null( K, scaleme )
% EMERGENCE_IO_NULL returns the marginal likelihood of a sequence under a
% fully random model that assume that
%   - "K": the length of the sequence to consider (note that we do not need
%       the actual sequence here).
%   - "scaleme": a string ('lin' or 'log') specifying whether the model
%       evidence sould be computed on a linear or logarithmic scale.
% 
% Copyright (c) 2018 Maxime Maheu

%% INITIALIZATION
%  ==============

% By default, return model evidence in log scale
if nargin < 2, scaleme = 'log'; end

%% SEQUENCE'S MARGINAL LIKELIHOOD
%  ==============================

% Compute the marginal likelihood of a sequence of length K under a 
% p(y) = (1/2)^N(A) + (1/2)^N(B)
% <=> p(y) = (1/2)^K
% <=> log(p(y)) = -K * log(2)
if     strcmpi(scaleme, 'lin'), pY = (1/2) .^ K;
elseif strcmpi(scaleme, 'log'), pY = -K .* log(2);
end

%% PREDICTION
%  ==========

% The likelihood that the next observation will be an A is simply chance
% level, which is in the case of binary sequences 1/2
if nargout > 1, pAgY = ones(1,numel(K)) ./ 2; end

end