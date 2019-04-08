function [Gxp] = cfp2per(Gx,p)
%cfp2per Convert cumulative frequency polygon to percentiles.
%   GXP = CFP2PER(GX,P) returns the percentiles of the cumulative frequency
%   polygon GX for the probabilities P. This function was adapted from the
%   code described in Appendix B, Ulrich et al. (2007).
%
%   See also RT2CFP, RT2CDF, RACEMODEL, SWITCHCOST, GETAUC.
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

% Get number of probabilities
np = length(p);

% Compute percentiles using linear interpolation
Gxp = zeros(np,1);
for i = 1:np
    idx = find(Gx<=p(i),1,'last');
    Gxp(i) = idx+(p(i)-Gx(idx))/(Gx(idx+1)-Gx(idx));
end