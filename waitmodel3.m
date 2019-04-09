function [Fx,Fy,Fz,Fxyz,Fwait,Fdiff,q] = waitmodel3(x,y,z,xyz,varargin)
%waitmodel3 Generate a wait model of trisensory reaction times.
%   [FX,FY,FZ,FXYZ] = WAITMODEL3(X,Y,Z,XYZ) returns the cumulative
%   distribution functions (CDFs) for the unisensory RT distributions X, Y
%   and Z, and the trisensory RT distribution XYZ at 10 linearly-spaced
%   quantiles. X, Y, Z and XYZ are not required to have an equal number of
%   observations. This function treats NaNs as missing values, and ignores
%   them.
%
%   To generate CDFs and wait models for the bisensory conditions XY, XZ
%   and YZ, use the function WAITMODEL on the corresponding unisensory and
%   bisensory RTs. To compare across bisensory and trisensory conditions,
%   use the same RT limits (see below).
%
%   [...,FWAIT] = WAITMODEL(...) returns the wait model based on the joint
%   probability of X, Y and Z (Townsend & Ashby, 1983). By default, the 
%   model assumes statistical independence between RTs on different sensory
%   channels, but this assumption can be specified using the DEP argument
%   (see below). For valid estimates of FWAIT, the stimuli used to generate
%   X, Y, Z and XYZ should be presented in random order to meet the
%   assumption of context invariance.
%
%   [...,FDIFF] = WAITMODEL3(...) returns the difference between FXYZ and
%   FWAIT to test for violations the wait model.
%
%   [...,Q] = WAITMODEL3(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test and the probabilities used to compute the
%   percentiles for the horizontal test.
%
%   [...] = WAITMODEL3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%               independence (AND model; default), -1 to assume perfect 
%               negative dependence (Diederich's bound) and 1 to assume
%               perfect positive dependence (Colonius's upper bound)
%   'test'      a string specifying how to test the wait model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%   'sharp'     a scalar specifying whether or not to sharpen the overly
%               conservative lower bound: pass in 1 to sharpen (Diederich's
%               bound; default) and 0 to not (Colonius's lower bound)
%
%   See also WAITMODEL, RACEMODEL3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Townsend JT, Eidels A (2011) Workload capacity spaces: A
%           unified methodology for response time measures of efficiency as
%           workload is varied. Psychon Bull Rev 18:659–681.
%       [2] Colonius H, Vorberg D (1994) Distribution inequalities for
%           parallel models with unlimited capacity. J Math Psychol
%           38:35-58.
%       [2] Diederich A (1992) Probability inequalities for testing
%           separate activation models of divided attention. Percept
%           Psychophys 14(2):247-279.
%       [3] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

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

% Compute wait model
if nargout > 4
    if dep == 0 % AND Model
        Fwait = Fx.*Fy.*Fz;
    elseif dep == -1
        if sharp == 1 % Diederich's Bound
            Fxy = Fx.*Fy; Fyz = Fy.*Fz;
            Fwait = min(Fxy+Fyz-Fy,ones(size(Fxyz)));
        elseif sharp == 0 % Colonius's Lower Bound
            Fwait = max(Fx+Fy+Fz-2,zeros(size(Fxyz)));
        end
    elseif dep == 1 % Colonius's Upper Bound
        Fwait = min([Fx,Fy,Fz],[],2);
    end
    if strcmpi(test,'hor')
        Fwait = cfp2per(Fwait,p);
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
        Fdiff = Fxyz-Fwait;
    elseif strcmpi(test,'hor')
        Fdiff = Fwait-Fxyz;
    end
end

% Get probabilities for horizontal test
if nargout > 5 &&  strcmpi(test,'hor')
    q = p;
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