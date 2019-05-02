function [Ccoef,Clim,Csup,q,lim] = andcapacity(x,y,xy,varargin)
%andcapacity Capacity coefficient for a bisensory AND task.
%   CCOEF = ANDCAPACITY(X,Y,XY) returns the capacity coefficient for a
%   bisensory AND task at 10 linearly-spaced quantiles. CCOEF values of 1
%   imply that the system has unlimited capacity, values below 1 imply
%   limited capacity and values above 1 imply super capacity (Townsend &
%   Eidels, 2011). X, Y and XY can have different lengths. This function
%   treats NaNs as missing values, and ignores them.
%
%   [...,CLIM,CSUP] = ANDCAPACITY(...) returns the predicted bounds of
%   extreme limited and super capacity, respectively.
%
%   [...,Q] = ANDCAPACITY(...) returns the RT quantiles used to compute the
%   CDFs.
%
%   [...,LIM] = ANDCAPACITY(...) returns the lower and upper RT limits used
%   to compute the CDFs. These values can be used to set the CDF limits of
%   subsequent tests that are to be compared with this one.
%
%   [...] = ANDCAPACITY(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'p'         a vector specifying the probabilities for computing the
%               quantiles of a vertical test or the percentiles of a
%               horizontal test (default=0.05:0.1:0.95)
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,XY]),max([X,Y,XY])])
%
%   See also ANDCAPACITY3, ORCAPACITY, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Townsend JT, Eidels A (2011) Workload capacity spaces: A
%           unified methodology for response time measures of efficiency as
%           workload is varied. Psychon Bull Rev 18:659–681.
%
%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 01-May-2019

% Decode input variable arguments
[p,lim] = decode_varargin(varargin);

% Transpose row vectors
if isrow(x), x = x'; end
if isrow(y), y = y'; end
if isrow(xy), xy = xy'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([x;y;xy]),max([x;y;xy])];
end

% Compute CDFs
Fx = rt2cdf(x,p,lim);
Fy = rt2cdf(y,p,lim);
[Fxy,q] = rt2cdf(xy,p,lim);

% Compute capacity coefficient
Ccoef = log(Fx.*Fy)./log(Fxy);

% Compute bounds of limited and super capacity
Clim = log(Fx.*Fy)./log(max(Fx+Fy-1,zeros(size(Fxy)))); % Colonius-Vorberg lower bound
Csup = log(Fx.*Fy)./log(min(Fx,Fy)); % Colonius-Vorberg upper bound

function [p,lim] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'lim')) && ~isempty(varargin{find(strcmpi(varargin,'lim'))+1})
    lim = varargin{find(strcmpi(varargin,'lim'))+1};
    if ~isnumeric(lim) || isscalar(lim) || any(isnan(lim)) || any(isinf(lim)) || any(lim<0) || lim(1)>=lim(2)
        error('LIM must be a 2-element vector of positive values.')
    end
else
    lim = []; % default: unspecified
end