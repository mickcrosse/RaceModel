function [Ccoef,Clim,Csup,q,lim] = andcapacity3(x,y,z,xyz,varargin)
%andcapacity3 Capacity coefficient for a trisensory AND task.
%   CCOEF = ANDCAPACITY3(X,Y,Z,XYZ) returns the capacity coefficient for a
%   trisensory AND task at 10 linearly-spaced quantiles.  CCOEF values of 1
%   imply that the system has unlimited capacity, values below 1 imply
%   limited capacity and values above 1 imply super capacity (Townsend &
%   Eidels, 2011). X, Y, Z and XYZ are not required to have an equal number
%   of observations. This function treats NaNs as missing values, and
%   ignores them.
%
%   [...,CLIM,CSUP] = ANDCAPACITY3(...) returns the predicted bounds of
%   limited and super capacity, respectively.
%
%   [...,Q] = ANDCAPACITY3(...) returns the RT quantiles used to compute
%   the CDFs.
%
%   [...,LIM] = ANDCAPACITY3(...) returns the lower and upper RT limits
%   used to compute the CDFs. These values can be used to set the CDF
%   limits of subsequent tests that are to be compared with this one.
%
%   [...] = ANDCAPACITY3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'p'         a vector specifying the probabilities for computing the
%               quantiles of a vertical test or the percentiles of a
%               horizontal test (default=0.05:0.1:0.95)
%   'outlier'   a 2-element vector specifying the lower and upper RT
%               cutoffs for outlier correction (default=no correction)
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'sharp'     a scalar specifying whether or not to sharpen the overly
%               conservative lower bound: pass in 1 to sharpen (Diederich's
%               bound; default) and 0 to not (Colonius-Vorberg lower bound)
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
%   Apr 2017; Last Revision: 15-Apr-2019

% Decode input variable arguments
[p,outlier,per,lim,sharp] = decode_varargin(varargin);

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
Fx = rt2cdf(x,p,lim);
Fy = rt2cdf(y,p,lim);
Fz = rt2cdf(z,p,lim);
[Fxyz,q] = rt2cdf(xyz,p,lim);

% Compute capacity coefficient
Ccoef = log(Fx.*Fy.*Fz)./log(Fxyz);

% Compute bounds of limited and super capacity
if sharp == 1
    Fxy = Fx.*Fy; Fyz = Fy.*Fz;
    Clim = log(Fx.*Fy.*Fz)./log(Fxy+Fyz-Fy); % Diederich's bound
elseif sharp == 0
    Clim = log(Fx.*Fy.*Fz)./abs(log(Fx+Fy+Fz-2)); % Colonius-Vorberg lower bound
end
Csup = log(Fx.*Fy.*Fz)./log(min([Fx,Fy,Fz],[],2)); % Colonius-Vorberg upper bound

function [p,outlier,per,lim,sharp] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'sharp')) && ~isempty(varargin{find(strcmpi(varargin,'sharp'))+1})
    sharp = varargin{find(strcmpi(varargin,'sharp'))+1};
    if sharp~=0 && sharp~=1
        error('SHARP must be a scalar with a value of 0 or 1.')
    end
else
    sharp = 1; % default: sharpen (Diederich's Bound)
end