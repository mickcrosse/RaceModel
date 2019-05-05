function [Ccoef,Clim,Csup,t] = andcapacity3(x,y,z,xyz,p,varargin)
%andcapacity3 Capacity coefficient for a trisensory AND task.
%   CCOEF = ANDCAPACITY3(X,Y,Z,XYZ) returns the capacity coefficient for a
%   trisensory AND task at 10 linearly-spaced intervals. CCOEF values of 1
%   imply that the system has unlimited capacity, values below 1 imply
%   limited capacity and values above 1 imply super capacity (Townsend &
%   Eidels, 2011). X, Y, Z and XYZ can have different lengths. This
%   function treats NaNs as missing values, and ignores them.
%
%   [...] = ANDCAPACITY3(...,P) uses the intervals P to compute CCOEF. P is
%   a vector of decimal values between 0 and 1 inclusive
%   (default=0.05:0.1:0.95).
%
%   [...,CLIM,CSUP] = ANDCAPACITY3(...) returns the predicted bounds of
%   extreme limited and super capacity, respectively.
%
%   [...,T] = ANDCAPACITY3(...) returns the time intervals used to compute
%   the CDFs for the vertical test.
%
%   [...] = ANDCAPACITY3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'sharp'     a scalar specifying whether or not to sharpen the overly
%               conservative lower bound: pass in 1 to sharpen (Diederich's
%               bound) and 0 to not (Colonius-Vorberg lower bound; default)
%
%   See also ANDCAPACITY, ORCAPACITY3, TPERMTEST, EFFECTSIZE.
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
[lim,sharp] = decode_varargin(varargin);

% Set default values
if ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
elseif nargin < 5 || isempty(p)
    p = 0.05:0.1:0.95;
end

% Transpose row vectors
if isrow(x), x = x'; end
if isrow(y), y = y'; end
if isrow(z), z = z'; end
if isrow(xyz), xyz = xyz'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([x;y;z;xyz]),max([x;y;z;xyz])];
end

% Compute CDFs
Fx = rt2cdf(x,p,lim);
Fy = rt2cdf(y,p,lim);
Fz = rt2cdf(z,p,lim);
[Fxyz,t] = rt2cdf(xyz,p,lim);

% Compute capacity coefficient
Ccoef = log(Fx.*Fy.*Fz)./log(Fxyz);

% Compute bounds of limited and super capacity
if sharp == 1
    Fxy = Fx.*Fy; Fyz = Fy.*Fz;
    Clim = log(Fx.*Fy.*Fz)./log(Fxy+Fyz-Fy); % Diederich's bound
elseif sharp == 0
    Clim = log(Fx.*Fy.*Fz)./log(max(Fx+Fy+Fz-2,zeros(size(Fxyz)))); % Colonius-Vorberg lower bound
end
Csup = log(Fx.*Fy.*Fz)./log(min([Fx,Fy,Fz],[],2)); % Colonius-Vorberg upper bound

function [lim,sharp] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'sharp')) && ~isempty(varargin{find(strcmpi(varargin,'sharp'))+1})
    sharp = varargin{find(strcmpi(varargin,'sharp'))+1};
    if sharp~=0 && sharp~=1
        error('SHARP must be a scalar with a value of 0 or 1.')
    end
else
    sharp = 0; % default: sharpen (Diederich's Bound)
end