function [gain] = rsegain(x,y,xy,varargin)
%rsegain Multisensory gain for a redundant signals effect.
%   GAIN = RSEGAIN(X,Y,XY) returns the multisensory gain for a redundant
%   signals effect (RSE), quantified by the area between the cumulative
%   distribution functions of the multisensory RT distribution XY and the
%   race model based on the unisensory RT distributions X and Y (Colonius
%   & Diederich, 2006). This function does not require X, Y and XY to have
%   the same number of observations. This function treats NaNs as missing
%   values, and ignores them.
%
%   [...] = RSEGAIN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'quant'     a vector specifying the quantiles used to compute the CDFs
%               (default=[0.05:0.05:1])
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'dep'       a scalar specifying whether statistical dependence between
%               X and Y is assumed: pass in 0 to assume independence (Raab,
%               1962; default), -1 to assume a perfect negative dependence
%               (Miller, 1982)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       entire area (default)
%                   'pos'       positive area
%                   'neg'       negative area
%
%   See also RACEMODEL, RSEBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Colonius H, Diederich A (2006) The Race Model Inequality:
%           Interpreting a Geometric Measure of the Amount of Violation.
%           Psychol Rev 113(1):148–154.
%       [2] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [3] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 6-Feb-2019

% ***Have not yet implemented horizontal test***

% Decode input variable arguments
[q,per,dep,test,area] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(3,2);
lims(1,:) = prctile(x,per);
lims(2,:) = prctile(y,per);
lims(3,:) = prctile(xy,per);

% Limit RTs to specified range
x = x(x>lims(1,1) & x<lims(1,2));
y = y(y>lims(2,1) & y<lims(2,2));
xy = xy(xy>lims(3,1) & xy<lims(3,2));

% Get min and max RT limits
lims = [min(lims(:)),max(lims(:))];

% Compute cumulative distribution functions
fx = rt2cdf(x,q,lims);
fy = rt2cdf(y,q,lims);
fxy = rt2cdf(xy,q,lims);

% Compute race model
if dep == 0 % Raab's Model
    frace = fx+fy-fx.*fy;
elseif dep == -1 % Miller's Bound
    frace = fx+fy;
    frace(frace>1) = 1;
end

% Compute difference
if strcmpi(test,'ver')
    fdiff = fxy-frace;
elseif strcmpi(test,'hor')
    fdiff = frace-fxy;
end

gain = getauc(q,fdiff,area);

function [auc] = getauc(x,y,p)
%getauc Get area under the curve.
%   Y = GETAUC(X,P) returns the area under the curve Y with respect to X
%   based on the portion of the curve P. Valid values for argument P are
%   'all' (entire portion), 'pos' (positive portion), and 'neg' (negative
%   portion).

% Replace negative/positive values with zeros
if strcmpi(p,'pos') % positive portion
    y(y<0) = 0;
elseif strcmpi(p,'neg') % negative portion
    y(y>0) = 0;
end

% Compute AUC using trapezoidal method
auc = trapz(x,y);

function [q,per,dep,test,area] = decode_varargin(varargin)
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
    end
else
    test = 'ver'; % default: vertical test
end
if any(strcmpi(varargin,'area')) && ~isempty(varargin{find(strcmpi(varargin,'area'))+1})
    area = varargin{find(strcmpi(varargin,'area'))+1};
    if ~any(strcmpi(area,{'all','pos'}))
        error('Invalid value for argument AREA. Valid values are: ''all'', ''pos''.')
    end
else
    area = 'all'; % default: entire area
end