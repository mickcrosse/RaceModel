function [Fx,Fy,Fz,Fxyz,Frace,Fdiff,q] = racemodel3(x,y,z,xyz,varargin)
%racemodel3 Generate a race model for trisensory reaction times.
%   [FX,FY,FZ,FXYZ] = RACEMODEL3(X,Y,Z,XYZ) returns the cumulative
%   distribution functions (CDFs) for the unisensory RT distributions X, Y
%   and Z, and the trisensory RT distribution XYZ at 10 linearly-spaced
%   quantiles. X, Y, Z and XYZ are not required to have an equal number of
%   observations. This function treats NaNs as missing values, and ignores
%   them.
%
%   To generate CDFs and race models for the bisensory conditions XY, XZ
%   and YZ, use the function RACEMODEL on the corresponding unisensory and
%   bisensory RTs. To compare across bisensory and trisensory conditions,
%   use the same RT limits (see below).
%
%   [...,FRACE] = RACEMODEL3(...) returns the race model based on
%   probability summation of X, Y and Z (Colonius et al., 2017). By
%   default, channel independence is assumed (Raab, 1962). For valid
%   estimates of FRACE, the stimuli used to generate X, Y, Z and XYZ should
%   be presented in random order to meet the assumption of context
%   invariance.
%
%   [...,FDIFF] = RACEMODEL3(...) returns the difference between FXYZ and
%   FRACE to test for violations the race model (Colonius et al., 2017).
%
%   [...,Q] = RACEMODEL3(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test.
%
%   [...] = RACEMODEL3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'p'         a vector specifying the probabilities for computing the
%               quantiles of a vertical test or the percentiles of a
%               horizontal test (default=0.05:0.1:0.95)
%   'outlier'   a 2-element vector specifying the lower and upper RT
%               cutoffs for outlier correction (default=no correction).
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly to other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'dep'       a scalar specifying the model's assumption of statistical
%               dependence between sensory channels: pass in 0 to assume
%               independence (Raab, 1962; default), -1 to assume a perfect
%               negative dependence (Miller, 1982) and 1 to assume a
%               perfect positive dependence (Grice et al., 1984)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%
%   See also RACEMODEL, RSEGAIN3, RSEBENEFIT3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Colonius H, Wolff FH, Diederich A (2017) Trimodal race model
%           inequalities in multisensory integration: I. Basics. Front
%           Psychol 8:1141.
%       [2] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [3] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.
%       [4] Grice GR, Canham L, Gwynne JW (1984) Absence of a redundant-
%           signals effect in a reaction time task with divided attention.
%           Percept Psychophys 36:565-570.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 4-Apr-2019

% Decode input variable arguments
[p,outlier,per,lim,dep,test] = decode_varargin(varargin);

% Outlier correction procedure
if ~isempty(outlier)
    x(x<outlier(1)|x>outlier(2)) = [];
    y(y<outlier(1)|y>outlier(2)) = [];
    z(z<outlier(1)|z>outlier(2)) = [];
    xyz(xyz<outlier(1)|xyz>outlier(2)) = [];
end

% Get RT range for each condition
lims = zeros(4,2);
lims(1,:) = prctile(x,per);
lims(2,:) = prctile(y,per);
lims(3,:) = prctile(z,per);
lims(4,:) = prctile(xyz,per);

% Limit RTs to specified range
x = x(x>=lims(1,1) & x<=lims(1,2));
y = y(y>=lims(2,1) & y<=lims(2,2));
z = z(z>=lims(3,1) & z<=lims(3,2));
xyz = xyz(xyz>=lims(4,1) & xyz<=lims(4,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
if strcmpi(test,'ver')
    Fx = rt2cdf(x,p,lim);
    Fy = rt2cdf(y,p,lim);
    Fz = rt2cdf(z,p,lim);
    [Fxyz,q] = rt2cdf(xyz,p,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fz = rt2cfp(z,lim(2));
    Fxyz = rt2cfp(xyz,lim(2));
end

% Compute race model
if nargout > 4
    if dep == 0 % Raab's Model
        Fxy = Fx+Fy-Fx.*Fy;
        Frace = Fxy+Fz-Fxy.*Fz;
    elseif dep == -1 % Miller's Bound
        Frace = min(Fx+Fy+Fz,ones(size(Fxyz)));
    elseif dep == 1 % Grice's Bound
        Frace = max([Fx,Fy,Fz],[],2);
    end
    if strcmpi(test,'hor')
        Frace = cfp2per(Frace,p);
    end
end

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    Fx = cfp2per(Fx,p);
    Fy = cfp2per(Fy,p);
    Fz = cfp2per(Fz,p);
    Fxyz = cfp2per(Fxyz,p);
end

% Compute difference
if nargout > 5
    if strcmpi(test,'ver')
        Fdiff = Fxyz-Frace;
    elseif strcmpi(test,'hor')
        Fdiff = Frace-Fxyz;
    end
end

function [p,outlier,per,lim,dep,test] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'outlier')) && ~isempty(varargin{find(strcmpi(varargin,'outlier'))+1})
    outlier = varargin{find(strcmpi(varargin,'outlier'))+1};
    if ~isnumeric(outlier) || isscalar(outlier) || any(isnan(outlier)) || any(isinf(outlier)) || any(outlier<0) || outlier(1)>=outlier(2)
        error('OUTLIER must be a 2-element vector of positive values.')
    end
else
    outlier = []; % default: unspecified
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