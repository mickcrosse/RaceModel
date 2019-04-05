function [Fx,Fy,Fxy,Frace,Fdiff,q] = racemodel(x,y,xy,varargin)
%racemodel Generate a race model of bisensory reaction times.
%   [FX,FY,FXY] = RACEMODEL(X,Y,XY) returns the cumulative distribution
%   functions (CDFs) for the unisensory RT distributions X and Y, and the
%   bisensory RT distribution XY at 10 linearly-spaced quantiles. X, Y and
%   XY are not required to have an equal number of observations. This
%   function treats NaNs as missing values, and ignores them.
%
%   [...,FRACE] = RACEMODEL(...) returns the race model based on
%   probability summation of X and Y. By default, the race model assumes
%   statistical independence between sensory channels (Raab, 1962). For
%   valid estimates of FRACE, the stimuli used to generate X, Y and XY
%   should be presented in random order to meet the assumption of context
%   invariance.
%
%   [...,FDIFF] = RACEMODEL(...) returns the difference between FXY and
%   FRACE to test for violations the race model (Miller, 1982).
%
%   [...,Q] = RACEMODEL(...) returns the quantiles used to compute the
%   CDFs.
%
%   [...] = RACEMODEL(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'p'         a vector specifying the probabilities for computing the
%               quantiles of a vertical test or the percentiles of a
%               horizontal test (default=0.05:0.1:0.95)
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly to other conditions
%               (default=[min([X,Y,XY]),max([X,Y,XY])])
%   'dep'       a scalar specifying the model's assumption of statistical
%               dependence between sensory channels: pass in 0 to assume
%               independence (Raab, 1962; default), -1 to assume a perfect
%               negative dependence (Miller, 1982) and 1 to assume a
%               perfect positive dependence (Grice et al., 1984)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%
%   See also RACEMODEL3, RSEGAIN, RSEBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [2] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.
%       [3] Grice GR, Canham L, Gwynne JW (1984) Absence of a redundant-
%           signals effect in a reaction time task with divided attention.
%           Percept Psychophys 36:565-570.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 4-Apr-2019

% Decode input variable arguments
[p,per,lim,dep,test] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(3,2);
lims(1,:) = prctile(x,per);
lims(2,:) = prctile(y,per);
lims(3,:) = prctile(xy,per);

% Limit RTs to specified range
x = x(x>=lims(1,1) & x<=lims(1,2));
y = y(y>=lims(2,1) & y<=lims(2,2));
xy = xy(xy>=lims(3,1) & xy<=lims(3,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
if strcmpi(test,'ver')
    Fx = rt2cdf(x,p,lim);
    Fy = rt2cdf(y,p,lim);
    [Fxy,q] = rt2cdf(xy,p,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fxy = rt2cfp(xy,lim(2));
end

% Compute race model
if nargout > 3
    if dep == 0 % Raab's Model
        Frace = Fx+Fy-Fx.*Fy;
    elseif dep == -1 % Miller's Bound
        Frace = min([Fx+Fy,ones(size(Fxy))],2);
    elseif dep == 1 % Grice's Bound
        Frace = max([Fx,Fy],[],2);
    end
    if strcmpi(test,'hor')
        Frace = cfp2per(Frace,p,lim(2));
    end
end

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    Fx = cfp2per(Fx,p,lim(2));
    Fy = cfp2per(Fy,p,lim(2));
    Fxy = cfp2per(Fxy,p,lim(2));
end

% Compute difference
if nargout > 4
    if strcmpi(test,'ver')
        Fdiff = Fxy-Frace;
    elseif strcmpi(test,'hor')
        Fdiff = Frace-Fxy;
    end
end

function [p,per,lim,dep,test] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'p')) && ~isempty(varargin{find(strcmpi(varargin,'p'))+1})
    p = varargin{find(strcmpi(varargin,'p'))+1};
    if ~isnumeric(p) || isscalar(p) || any(isnan(p)) || any(isinf(p)) || any(p<0) || any(p>1) || any(diff(p)<=0)
        error('P must be a vector with values between 0 and 1.')
    end
else
    p = 0.05:0.1:0.95; % default: 0.05 to 0.95 in 0.1 increments
end
if any(strcmpi(varargin,'per')) && ~isempty(varargin{find(strcmpi(varargin,'per'))+1})
    per = varargin{find(strcmpi(varargin,'per'))+1};
    if ~isnumeric(per) || isscalar(per) || any(isnan(per)) || any(isinf(per)) || any(per<0) || any(per>100) || per(1)>=per(2)
        error('PER must be a 2-element vector with values between 0 and 100.')
    end
else
    per = [0,100]; % default: all RTs
end
if any(strcmpi(varargin,'lim')) && ~isempty(varargin{find(strcmpi(varargin,'lim'))+1})
    lim = varargin{find(strcmpi(varargin,'lim'))+1};
    if ~isnumeric(lim) || isscalar(lim) || any(isnan(lim)) || any(isinf(lim)) || any(lim<0) || lim(1)>=lim(2)
        error('LIM must be a 2-element vector of positive values.')
    end
else
    lim = []; % default: unspecified
end
if any(strcmpi(varargin,'dep')) && ~isempty(varargin{find(strcmpi(varargin,'dep'))+1})
    dep = varargin{find(strcmpi(varargin,'dep'))+1};
    if dep~=0 && dep~=-1 && dep~=1
        error('DEP must be a scalar with a value of -1, 0 or 1.')
    end
else
    dep = 0; % default: assume statistical independence (Raab's Model)
end
if any(strcmpi(varargin,'test')) && ~isempty(varargin{find(strcmpi(varargin,'test'))+1})
    test = varargin{find(strcmpi(varargin,'test'))+1};
    if ~any(strcmpi(test,{'ver','hor'}))
        error('Invalid value for argument TEST. Valid values are: ''ver'', ''hor''.')
    end
else
    test = 'ver'; % default: vertical test
end