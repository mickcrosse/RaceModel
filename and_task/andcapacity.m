function [Ccoef,Clim,Csup,t] = andcapacity(x,y,xy,p,varargin)
%andcapacity Capacity coefficient for a bisensory AND task.
%   CCOEF = ANDCAPACITY(X,Y,XY) returns the capacity coefficient for a
%   bisensory AND task at 10 linearly-spaced intervals. CCOEF values of 1
%   imply that the system has unlimited capacity, values below 1 imply
%   limited capacity and values above 1 imply super capacity (Townsend &
%   Eidels, 2011). X, Y and XY can have different lengths. This function
%   treats NaNs as missing values, and ignores them.
%
%   [...] = ANDCAPACITY(...,P) uses the intervals P to compute CCOEF. P is
%   a vector of decimal values between 0 and 1 inclusive
%   (default=0.05:0.1:0.95).
%
%   [...,CLIM,CSUP] = ANDCAPACITY(...) returns the predicted bounds of
%   extreme limited and super capacity, respectively.
%
%   [...,T] = ANDCAPACITY(...) returns the time intervals used to compute
%   the CDFs for the vertical test.
%
%   [...] = ANDCAPACITY(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
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
%   Apr 2017; Last Revision: 3-May-2019

% Decode input variable arguments
lim = decode_varargin(varargin);

% Set default values
if nargin < 4 || isempty(p)
    p = 0.05:0.1:0.95;
elseif ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
end

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
[Fxy,t] = rt2cdf(xy,p,lim);

% Compute capacity coefficient
Ccoef = log(Fx.*Fy)./log(Fxy);

% Compute bounds of limited and super capacity
Clim = log(Fx.*Fy)./log(max(Fx+Fy-1,zeros(size(Fxy)))); % Colonius-Vorberg lower bound
Csup = log(Fx.*Fy)./log(min(Fx,Fy)); % Colonius-Vorberg upper bound

function lim = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'lim')) && ~isempty(varargin{find(strcmpi(varargin,'lim'))+1})
    lim = varargin{find(strcmpi(varargin,'lim'))+1};
    if ~isnumeric(lim) || isscalar(lim) || any(isnan(lim)) || any(isinf(lim)) || any(lim<0) || lim(1)>=lim(2)
        error('LIM must be a 2-element vector of positive values.')
    end
else
    lim = []; % default: unspecified
end