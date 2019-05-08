function [Fx,Fy,Fz,Fxyz,Fmodel,t] = ormodel3(x,y,z,xyz,p,varargin)
%ormodel3 Model trisensory reaction times for an OR task.
%   [FX,FY,FZ,FXYZ] = ORMODEL3(X,Y,Z,XYZ) returns the cumulative
%   distribution functions (CDFs) for the unisensory RT distributions X, Y
%   and Z, and the trisensory RT distribution XYZ at 10 linearly-spaced
%   intervals. X, Y, Z and XYZ can have different lengths. This function
%   treats NaNs as missing values, and ignores them.
%
%   To compute CDFs and models for the bisensory conditions XY, XZ and YZ,
%   use the function ORMODEL on the corresponding unisensory and bisensory
%   RTs. To compare across bisensory and trisensory conditions, use the
%   same RT limits (see below).
%
%   [...] = ORMODEL3(...,P) uses the intervals P to generate CDFs. P is a
%   vector of decimal values between 0 and 1 inclusive. For horizontal
%   tests, P is the probabilities used to compute the CDF quantiles
%   (default=0.05:0.1:0.95).
%
%   [...,FMODEL] = ORMODEL3(...) returns the OR (race) model based on the
%   probability summation of X, Y and Z (Raab, 1962). By default, the model
%   assumes stochastic independence between RTs on different sensory
%   channels, but this assumption can be specified using the DEP argument
%   (see below). For valid estimates of FMODEL, the stimuli used to
%   generate X, Y, Z and XYZ should be presented in random order to meet
%   the assumption of context invariance.
%
%   [...,T] = ORMODEL3(...) returns the time intervals used to compute the
%   CDFs for the vertical test.
%
%   [...] = ORMODEL3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%               independence (OR model; default), -1 to assume perfect
%               negative dependence (Diederich's bound) and 1 to assume
%               perfect positive dependence (Grice's bound)
%   'test'      a string specifying how to test the model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%   'sharp'     a scalar specifying whether or not to sharpen the overly
%               conservative upper bound: pass in 1 to sharpen (Diederich's
%               bound; default) and 0 to not (Miller's bound)
%
%   See also ORMODEL, ORGAIN3, ORBENEFIT3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [3] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.
%       [4] Diederich A (1992) Probability inequalities for testing
%           separate activation models of divided attention. Percept
%           Psychophys 14(2):247-279.
%       [5] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 3-May-2019

% Decode input variable arguments
[lim,dep,test,sharp] = decode_varargin(varargin);

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
if nargout > 4
    if dep == 0 % OR model
        Fxy = Fx+Fy-Fx.*Fy;
        Fmodel = Fxy+Fz-Fxy.*Fz;
    elseif dep == -1
        Fupr = ones(length(p),1);
        if sharp == 1 % Diederich's bound
            Fxy = Fx+Fy-Fx.*Fy; Fxz = Fx+Fz-Fx.*Fz; Fyz = Fy+Fz-Fy.*Fz;
            F1 = Fxy+Fxz-Fx; F2 = Fxy+Fyz-Fy; F3 = Fxz+Fyz-Fz;
            Fmodel = min([F1,F2,F3,Fupr],[],2);
        elseif sharp == 0 % Miller's bound
            Fmodel = min(Fx+Fy+Fz,Fupr);
        end
    elseif dep == 1 % Grice's bound
        Fmodel = max([Fx,Fy,Fz],[],2);
    end
    if strcmpi(test,'hor')
        Fmodel = cfp2q(Fmodel,p);
    end
end

% Compute quantiles for horizontal test
if strcmpi(test,'hor')
    Fx = cfp2q(Fx,p);
    Fy = cfp2q(Fy,p);
    Fz = cfp2q(Fz,p);
    Fxyz = cfp2q(Fxyz,p);
end

% Time intervals for horizontal test not required
if nargout > 5 && strcmpi(test,'hor')
    error('Time intervals T not required for horizontal test.')
end

function [lim,dep,test,sharp] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'sharp')) && ~isempty(varargin{find(strcmpi(varargin,'sharp'))+1})
    sharp = varargin{find(strcmpi(varargin,'sharp'))+1};
    if sharp~=0 && sharp~=1
        error('SHARP must be a scalar with a value of 0 or 1.')
    end
else
    sharp = 1; % default: sharpen (Diederich's Bound)
end