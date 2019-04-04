function [Fx] = cfp2per(Gx,q,rtmax)
%cfp2per Convert cumulative frequency polygon to percentiles.
%   FX = CFP2PER(GX,Q,RTMAX) returns the percentiles of the cumulative
%   frequency polygon GX for the quantiles Q between the minimum RT value
%   and RTMAX as described in Ulrich et al. (2007). This funciton was
%   adapted from the code described in Appendix B, Ulrich et al. (2007).
%
%   See also RT2CFP, RT2CDF, RACEMODEL, SWITCHCOST, GETAUC.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res Meth
%           39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 3-Apr-2019

% Round max RT value
rtmax = round(rtmax);

% Compute RT value at each percentile
Fx = zeros(length(q),1);
for i = 1:length(q)
    clim = 100;
    for t = 1:rtmax
        if abs(Gx(t)-q(i)) < clim
            c = t;
            clim = abs(Gx(t)-q(i));
        end
    end
    if q(i) > Gx(c)
        Fx(i) = c+(q(i)-Gx(c))/(Gx(c+1)-Gx(c));
    else
        Fx(i) = c+(q(i)-Gx(c))/(Gx(c)-Gx(c-1));
    end
end