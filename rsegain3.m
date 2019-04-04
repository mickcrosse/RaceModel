function gain = rsegain3(x,y,z,xyz,varargin)
%rsegain3 Multisensory gain of a trisensory redundant signals effect.
%   GAIN = RSEGAIN3(X,Y,Z,XYZ) returns the multisensory gain of a redundant
%   signals effect (RSE), quantified by the area between the cumulative
%   distribution functions of the trisensory RT distribution XYZ, and the
%   race model based on the unisensory RT distributions X, Y and Z (Miller,
%   1986; Colonius & Diederich, 2006). This function does not require X, Y,
%   Z and XYZ to have an equal number of observations. This function treats
%   NaNs as missing values, and ignores them.
%
%   To compute the gain for the bisensory conditions XY, XZ and YZ, use the
%   function RSEGAIN on the corresponding unisensory and bisensory RTs. To
%   compare across bisensory and trisensory conditions, use the same lower
%   and upper RT limits (see below).
%
%   [...] = RSEGAIN3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'q'         a vector specifying the quantiles to be used to compute the
%               CDFs (default=[0.05:0.05:1])
%   'per'       a 2-element vector specifying the lower and upper RT
%               percentiles to be used for each condition (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               to be used to compute the CDFs: it is recommended to leave
%               this unspecified unless comparing directly to other
%               conditions (default=[min([x,y,z,xyz]),max([x,y,z,xyz])])
%   'dep'       a scalar specifying the model's assumption of statistical
%               dependence between X and Y: pass in 0 to assume
%               independence (Raab, 1962; default) or -1 to assume a
%               perfect negative dependence (Miller, 1982)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%
%   See also RSEGAIN, RACEMODEL3, RSEBENEFIT3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Miller J (1986) Timecourse of coactivation in bimodal divided
%           attention. Percept Psychophys 40(5):331-343.
%       [2] Colonius H, Diederich A (2006) The Race Model Inequality:
%           Interpreting a Geometric Measure of the Amount of Violation.
%           Psychol Rev 113(1):148–154.
%       [3] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [4] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 2-Apr-2019

% Decode input variable arguments
[q,per,lim,dep,test,area] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(4,2);
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

% Compute CDFs
if strcmpi(test,'ver')
    Fx = rt2cdf(x,q,lim);
    Fy = rt2cdf(y,q,lim);
    Fz = rt2cdf(z,q,lim);
    Fxyz = rt2cdf(xyz,q,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fz = rt2cfp(z,lim(2));
    Fxyz = rt2cfp(xyz,lim(2));
end

% Compute race model
if dep == 0 % Raab's Model
    Fxy = Fx+Fy-Fx.*Fy;
    Frace = Fxy+Fz-Fxy.*Fz;
elseif dep == -1 % Miller's Bound
    Frace = Fx+Fy+Fz;
end

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    Fxyz = cfp2per(Fxyz,q,lim(2));
    Frace = cfp2per(Frace,q,lim(2));
end

% Normalize race model between 0 and 1
Frace(Frace>1) = 1;

% Compute difference
if strcmpi(test,'ver')
    Fdiff = Fxyz-Frace;
elseif strcmpi(test,'hor')
    Fdiff = Frace-Fxyz;
end

% Compute multisensory gain
gain = getauc(q,Fdiff,area);

function [q,per,lim,dep,test,area] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'q')) && ~isempty(varargin{find(strcmpi(varargin,'q'))+1})
    q = varargin{find(strcmpi(varargin,'q'))+1};
    if ~isnumeric(q) || isscalar(q) || any(isnan(q)) || any(isinf(q)) || any(q<0) || any(q>1) || any(diff(q)<=0)
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
    elseif strcmpi(test,'hor')
        error('Horizontal test not yet implemented. Please watch out for updates.')
    end
else
    test = 'ver'; % default: vertical test
end
if any(strcmpi(varargin,'area')) && ~isempty(varargin{find(strcmpi(varargin,'area'))+1})
    area = varargin{find(strcmpi(varargin,'area'))+1};
    if ~any(strcmpi(area,{'all','pos','neg'}))
        error('Invalid value for argument AREA. Valid values are: ''all'', ''pos'', ''neg''.')
    end
else
    area = 'all'; % default: use all values
end