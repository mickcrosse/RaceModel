function [gain,Fdiff,t] = andgain3(x,y,z,xyz,p,varargin)
%andgain3 Multisensory gain for a trisensory AND task.
%   GAIN = ANDGAIN3(X,Y,Z,XYZ) returns the multisensory gain for a
%   trisensory AND task, quantified by the area between the CDFs of the
%   trisensory RT distribution XYZ, and the AND model based on the joint
%   probability of  X, Y and Z (Crosse et al., 2019). X, Y, Z and XYZ can
%   have different lengths. This function treats NaNs as missing values,
%   and ignores them.
%
%   To compute multisensory gain for the bisensory conditions XY, XZ and
%   YZ, use the function ANDGAIN on the corresponding unisensory and
%   bisensory RTs. To compare across bisensory and trisensory conditions,
%   use the same RT limits (see below).
%
%   [...] = ANDGAIN3(...,P) uses the intervals P to generate CDFs. P is a
%   vector of decimal values between 0 and 1 inclusive. For horizontal
%   tests, P is the probabilities used to compute the CDF quantiles
%   (default=0.05:0.1:0.95).
%
%   [...,FDIFF] = ANDGAIN3(...) returns the difference at each quantile to
%   test for violations of the model.
%
%   [...,T] = ANDGAIN3(...) returns the time intervals used to compute the
%   CDFs for the vertical test.
%
%   [...] = ANDGAIN3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'dep'       a scalar specifying the model's assumption of stochastic
%               dependence between sensory channels: pass in 0 to assume
%               independence (AND model; default), -1 to assume perfect
%               negative dependence (Diederich's bound) and 1 to assume
%               perfect positive dependence (Colonius-Vorberg upper bound)
%   'test'      a string specifying how to test the model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%   'sharp'     a scalar specifying whether or not to sharpen the overly
%               conservative upper bound: pass in 1 to sharpen (Diederich's
%               bound; default) and 0 to not (Colonius-Vorberg lower bound)
%
%   See also ANDGAIN, ANDMODEL3, ANDBENEFIT3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Colonius H, Diederich A (2006) The Race Model Inequality:
%           Interpreting a Geometric Measure of the Amount of Violation.
%           Psychol Rev 113(1):148–154.
%       [3] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 3-May-2019

% Decode input variable arguments
[lim,dep,test,area,sharp] = decode_varargin(varargin);

% Set default values
if nargin < 5 || isempty(p)
    p = 0.05:0.1:0.95;
elseif ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
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
if strcmpi(test,'ver')
    Fx = rt2cdf(x,p,lim);
    Fy = rt2cdf(y,p,lim);
    Fz = rt2cdf(z,p,lim);
    [Fxyz,t] = rt2cdf(xyz,p,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fz = rt2cfp(z,lim(2));
    Fxyz = rt2cfp(xyz,lim(2));
end

% Compute model
if dep == 0 % AND model
    Fmodel = Fx.*Fy.*Fz;
elseif dep == -1
    Flwr = zeros(length(p),1);
    if sharp == 1 % Diederich's bound
        Fxy = Fx.*Fy; Fxz = Fx.*Fz; Fyz = Fy.*Fz;
        F1 = Fxy+Fxz-Fx; F2 = Fxy+Fyz-Fy; F3 = Fxz+Fyz-Fz;
        Fmodel = max([F1,F2,F3,Flwr],[],2);
    elseif sharp == 0 % Colonius-Vorberg lower bound
        Fmodel = max(Fx+Fy+Fz-2,Flwr);
    end
elseif dep == 1 % Colonius-Vorberg upper bound
    Fmodel = min([Fx,Fy,Fz],[],2);
end

% Compute quantiles for horizontal test
if strcmpi(test,'hor')
    Fxyz = cfp2q(Fxyz,p);
    Fmodel = cfp2q(Fmodel,p);
end

% Compute difference
if strcmpi(test,'ver')
    Fdiff = Fxyz-Fmodel;
elseif strcmpi(test,'hor')
    Fdiff = Fmodel-Fxyz;
end

% Compute multisensory gain
gain = getauc(p,Fdiff,area);

% Time intervals for horizontal test not required
if nargout > 2 && strcmpi(test,'hor')
    error('Time intervals T not required for horizontal test.')
end

function [lim,dep,test,area,sharp] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'dep')) && ~isempty(varargin{find(strcmpi(varargin,'dep'))+1})
    dep = varargin{find(strcmpi(varargin,'dep'))+1};
    if dep~=-1 && dep~=0 && dep~=1
        error('DEP must be a scalar with a value of -1, 0 or 1.')
    end
else
    dep = 0; % default: assume stochastic independence
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
    if ~any(strcmpi(area,{'all','pos','neg'}))
        error('Invalid value for argument AREA. Valid values are: ''all'', ''pos'', ''neg''.')
    end
else
    area = 'all'; % default: use all values
end
if any(strcmpi(varargin,'sharp')) && ~isempty(varargin{find(strcmpi(varargin,'sharp'))+1})
    sharp = varargin{find(strcmpi(varargin,'sharp'))+1};
    if sharp~=0 && sharp~=1
        error('SHARP must be a scalar with a value of 0 or 1.')
    end
else
    sharp = 1; % default: sharpen (Diederich's Bound)
end