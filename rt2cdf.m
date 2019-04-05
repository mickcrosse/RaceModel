function [Fx,q] = rt2cdf(x,p,lim)
%rt2cdf Convert reaction times to cumulative probabilities.
%   FX = RT2CDF(X,P,LIM) returns the cumulative distribution function (CDF)
%   of the RT distribution X for the set of quantiles between the RT limits 
%   LIM corresponding to the probabilities P. P must be a vector of values 
%   between 0 and 1 and LIM must be a 2-element vector containing the lower 
%   and upper RT limits. CDFs drawn from different samples can be averaged 
%   as long as the same probabilities were used to compute the quantiles 
%   (Ratcliff, 1979). This function treats NaNs as missing values, and 
%   ignores them.
% 
%   [...,Q] = RT2CDF(...) returns the quantiles used to compute the CDF.
%
%   See also RT2CFP, CFP2PER, RACEMODEL, SWITCHCOST, GETAUC.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Ratcliff R (1979) Group reaction time distributions and an
%           analysis of distribution statistics. Psychol Bull 86(3):446-461.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 4-Apr-2019

% Set default values
if nargin < 3 || isempty(lim)
    lim = [min(x),max(x)];
elseif lim(1) > lim(2)
    error('Value of LIM(1) must be < LIM(2).')
end
if nargin < 2 || isempty(p)
    p = 0.05:0.1:0.95;
end

% Get number of observations
nx = sum(~isnan(x));
nq = length(p);

% Compute quantiles
q = lim(1)+p*(lim(2)-lim(1));

% Compute CDF
Fx = zeros(nq,1);
for i = 1:nq
    Fx(i) = sum(x<=q(i))/nx;
end