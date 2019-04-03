function fx = rt2cdf(x,q,lim)
%rt2cdf Convert reaction times to cumulative probabilities.
%   FX = RT2CDF(X,Q,LIM) returns the cumulative distribution function (CDF)
%   of the RT distribution X for quantiles Q between the RT limits LIM. Q
%   must be a vector of linearly-spaced quantiles between 0 and 1. LIM must
%   be a 2-element vector containing the lower and upper RT limits of the
%   CDF. CDFs drawn from different samples can be averaged as long as the
%   values of Q are kept constant (Ratcliff, 1979). This function treats
%   NaNs as missing values, and ignores them.
%
%   See also RACEMODEL, SWITCHCOST, GETAUC.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Ratcliff R (1979) Group reaction time distributions and an
%           analysis of distribution statistics. Psychol Bull 86(3):446-461.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 6-Feb-2019

% Set default values
if nargin < 3 || isempty(lim)
    lim = [min(x),max(x)];
elseif lim(1) > lim(2)
    error('Value of LIM(1) must be < LIM(2).')
end
if nargin < 2 || isempty(q)
    q = 0.05:0.05:1;
end

% Get number of observations
nx = sum(~isnan(x));
nq = length(q);

% Compute linearly-spaced quantiles between RT limits
qntls = linspace(lim(1),lim(2),nq);

% Compute CDF
fx = zeros(nq,1);
for i = 1:length(q)
    fx(i) = sum(x<=qntls(i))/nx;
end