function Gx = rt2cfp(x,tmax)
%rt2cfp Convert reaction times to a cumulative frequency polygon.
%   GX = RT2CFP(X,RTMAX) returns the cumulative frequency polygon of the RT
%   distribution X for every integer value between zero and TMAX. GX is
%   obtained by generating a step function using the sorted values of X,
%   calculating the midpoint of each vertical step segment and linearly
%   interpolating between adjacent midpoints. For a mathematical
%   description, see Appendix A, Ulrich et al. (2007).
%
%   See also CFP2PER, RT2CDF, RT2PDF, GETAUC.
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
%   Apr 2017; Last Revision: 3-May-2019

% Sort and Round RTs
x = round(sort(x));
tmax = round(tmax);

% Identify unique observations
[~,idx] = unique(x);

% Get number of replications
n = hist(x,x(idx));

% Compute cumulative sum
s = cumsum(n);

% Get rid of ties
x = x(idx);

% Get number of unique observations
k = length(x);

% Determine midpoints of vertical step segments
Gx = zeros(tmax,1);
for i = 1:k-1
    for t = x(i):x(i+1)-1
        if i == 1
            Gx(t) = 1/s(k) * (n(i)/2 + (n(i)+n(i+1))/2 * (t-x(i))/(x(i+1)-x(i)));
        else
            Gx(t) = 1/s(k) * (s(i-1) + n(i)/2 + (n(i)+n(i+1))/2 * (t-x(i))/(x(i+1)-x(i)));
        end
    end
end
Gx(x(k):tmax) = 1;