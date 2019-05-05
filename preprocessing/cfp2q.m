function q = cfp2q(Gx,p)
%cfp2q Convert cumulative frequency polygon to quantiles.
%   Q = CFP2Q(GX,P) returns the quantiles of the cumulative frequency
%   polygon GX for the probabilities P. This function was adapted from the
%   code described in Appendix B, Ulrich et al. (2007).
%
%   See also RT2CFP, RT2CDF, RT2PDF, GETAUC.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 4-Apr-2019

% Set default values
if nargin < 2 || isempty(p)
    p = 0.05:0.1:0.95;
end

% Get number of quantiles
nq = length(p);

% Compute quantiles using linear interpolation
q = zeros(nq,1);
for i = 1:nq
    idx = find(Gx<=p(i),1,'last');
    q(i) = idx+(p(i)-Gx(idx))/(Gx(idx+1)-Gx(idx));
end