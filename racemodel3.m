function [fx,fy,fz,fxyz,frace,fdiff] = racemodel3(x,y,z,xyz,varargin)
%racemodel3 Generate a race model based on trimodal reaction times.
%   [FX,FY,FZ,FXYZ] = RACEMODEL3(X,Y,Z,XYZ) returns the cumulative
%   distribution functions (CDFs) for the unisensory RT distributions X, Y
%   and Z, and the multisensory RT distribution XYZ at 20 linearly-spaced
%   quantiles between 0.05 and 1. This function does not require X, Y, Z
%   and XYZ to have an equal number of observations. This function treats
%   NaNs as missing values, and ignores them.
%
%   To generate CDFs and race models for conditions XY, XZ and YZ, use the
%   function RACEMODEL separately for each combination of unisensory and
%   multisensory RTs and input the same lower and upper values for argument
%   LIM to compare across datasets (see below).
%
%   [...,FRACE] = RACEMODEL3(...) returns the CDF of the race model based
%   on the unisensory RT distributions X, Y and Z (Colonius et al., 2017).
%   The race model is computed using probability summation (Raab, 1962),
%   which assumes statistical independence between X, Y and Z. For valid
%   estimates of FRACE, the stimuli used to generate X, Y, Z and XYZ should
%   be randomly interleaved in order to uphold the assumption of context
%   invariance (Luce, 1986).
%
%   [...,FDIFF] = RACEMODEL3(...) returns the difference between FXYZ and
%   FRACE to test whether XYZ exceeded statistical facilitation predicted
%   by the race model (Miller, 1982).
%
%   [...] = RACEMODEL3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'q'         a vector specifying the quantiles used to compute the CDFs
%               (default=[0.05:0.05:1])
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               used to compute the CDFs: it is recommended to leave this
%               unspecified or empty unless comparing to other conditions
%   'dep'       a scalar specifying whether statistical dependence between
%               X, Y and Z is assumed: pass in 0 to assume independence
%               (Raab, 1962; default), -1 to assume a perfect negative
%               dependence (Miller, 1982) and 1 to assume a perfect
%               positive dependence (Grice et al., 1986)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%
%   See also RACEMODEL, RSEGAIN, RSEBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Colonius H, Wolff FH, Diederich A (2017) Trimodal race model
%           inequalities in multisensory integration: I. Basics. Front
%           Psychol 8:1141.
%       [2] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [3] Luce RD (1986) Response times: Their role in inferring mental
%           organization. New York, NY: Oxford University Press.
%       [4] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.
%       [5] Grice GR, Canham L, Gwynne JW (1984) Absence of a redundant-
%           signals effect in a reaction time task with divided attention.
%           Percept Psychophys 36:565-570.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 7-Feb-2019

% Decode input variable arguments
[q,per,lim,dep,test] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(3,2);
lims(1,:) = prctile(x,per);
lims(2,:) = prctile(y,per);
lims(3,:) = prctile(z,per);
lims(4,:) = prctile(xyz,per);

% Limit RTs to specified range
x = x(x>lims(1,1) & x<lims(1,2));
y = y(y>lims(2,1) & y<lims(2,2));
z = z(z>lims(3,1) & z<lims(3,2));
xyz = xyz(xyz>lims(4,1) & xyz<lims(4,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute cumulative distribution functions
fx = rt2cdf(x,q,lim);
fy = rt2cdf(y,q,lim);
fz = rt2cdf(z,q,lim);
fxyz = rt2cdf(xyz,q,lim);

% Compute race model
if nargout > 3
    if dep == 0 % Raab's Model
        fxy = fx+fy-fx.*fy;
        frace = fxy+fz-fxy.*fz;
    elseif dep == -1 % Miller's Bound
        frace = fx+fy+fz;
        frace(frace>1) = 1;
    elseif dep == 1 % Grice's Bound
        frace = max([fx,fy,fz],[],2);
    end
end

% Compute difference
if nargout > 4
    if strcmpi(test,'ver')
        fdiff = fxyz-frace;
    elseif strcmpi(test,'hor')
        fdiff = frace-fxyz;
    end
end

function [q,per,lim,dep,test] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'q')) && ~isempty(varargin{find(strcmpi(varargin,'q'))+1})
    q = varargin{find(strcmpi(varargin,'q'))+1};
    if ~isnumeric(q) || isscalar(q) || any(isnan(q)) || any(isinf(q)) || any(q<0) || any(q>1) || q(1)>=q(2)
        error('Q must be a vector with values between 0 and 1.')
    end
else
    q = 0.05:0.05:1; % default: 0.05 to 1 in 0.05 increments
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
    if any(dep~=0) && any(dep~=-1) && any(dep~=1)
        error('DEP must be a scalar with a value of -1, 0 or 1.')
    end
else
    dep = 0; % default: assume statistical independence (Raab's Model)
end
if any(strcmpi(varargin,'test')) && ~isempty(varargin{find(strcmpi(varargin,'test'))+1})
    test = varargin{find(strcmpi(varargin,'test'))+1};
    if ~any(strcmpi(test,{'ver','hor'}))
        error('Invalid value for argument TEST. Valid values are: ''ver'', ''hor''.')
    elseif strcmpi(test,'hor')
        error('Horizontal test not yet implemented. Please watch out for updates.')
    end
else
    test = 'ver'; % default: vertical test
end