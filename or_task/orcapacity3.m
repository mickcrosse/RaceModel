function [Ccoef,Clim,Csup,q,lim] = orcapacity3(x,y,z,xyz,varargin)
%orcapacity3 Capacity coefficient for a trisensory OR task.
%   CCOEF = ORCAPACITY3(X,Y,Z,XYZ) returns the capacity coefficient for a
%   trisensory OR task at 10 linearly-spaced quantiles. CCOEF values of 1
%   imply that the system has unlimited capacity, values below 1 imply
%   limited capacity and values above 1 imply super capacity (Townsend &
%   Eidels, 2011). X, Y, Z and XYZ can have different lengths. This
%   function treats NaNs as missing values, and ignores them.
%
%   [...,CLIM,CSUP] = ORCAPACITY3(...) returns the predicted bounds of
%   extreme limited and super capacity, respectively.
%
%   [...,Q] = ORCAPACITY3(...) returns the RT quantiles used to compute the
%   CDFs.
%
%   [...,LIM] = ORCAPACITY3(...) returns the lower and upper RT limits used
%   to compute the CDFs. These values can be used to set the CDF limits of
%   subsequent tests that are to be compared with this one.
%
%   [...] = ORCAPACITY3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'sharp'     a scalar specifying whether or not to sharpen the overly
%               conservative upper bound: pass in 1 to sharpen (Diederich's
%               bound) and 0 to not (Miller's bound; default)
%
%   See also ORCAPACITY, ANDCAPACITY3, TPERMTEST, EFFECTSIZE.
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
[p,lim,sharp] = decode_varargin(varargin);

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
[Fxyz,q] = rt2cdf(xyz,p,lim);

% Compute survivor functions
Sx = 1-Fx;
Sy = 1-Fy;
Sz = 1-Fz;
Sxyz = 1-Fxyz;

% Compute capacity coefficient
Ccoef = log(Sxyz)./log(Sx.*Sy.*Sz);

% Compute bounds of limited and super capacity
Clim = log(min([Sx,Sy,Sz],[],2))./log(Sx.*Sy.*Sz); % Grice's bound
if sharp == 1
    Sxy = Sx+Sy-Sx.*Sy; Syz = Sy+Sz-Sy.*Sz;
    Csup = log(Sxy+Syz-Sy)./log(Sx.*Sy.*Sz); % Diederich's bound
elseif sharp == 0
    Csup = log(max(Sx+Sy+Sz-2,zeros(size(Sxyz))))./log(Sx.*Sy.*Sz); % Miller's bound
end

function [p,lim,sharp] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'sharp')) && ~isempty(varargin{find(strcmpi(varargin,'sharp'))+1})
    sharp = varargin{find(strcmpi(varargin,'sharp'))+1};
    if sharp~=0 && sharp~=1
        error('SHARP must be a scalar with a value of 0 or 1.')
    end
else
    sharp = 0; % default: sharpen (Diederich's Bound)
end