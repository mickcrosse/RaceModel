function [Fx,Fy,Fxy,Fmodel,t] = andmodel(x,y,xy,p,varargin)
%andmodel Model bisensory reaction times for an AND task.
%   [FX,FY,FXY] = ANDMODEL(X,Y,XY) returns the cumulative distribution
%   functions (CDFs) for the unisensory RT distributions X and Y, and the
%   bisensory RT distribution XY at 10 linearly-spaced intervals. X, Y and
%   XY can have different lengths. This function treats NaNs as missing
%   values, and ignores them.
%
%   [...] = ANDMODEL(...,P) uses the intervals P to generate CDFs. P is a
%   vector of decimal values between 0 and 1 inclusive. For horizontal
%   tests, P is the probabilities used to compute the CDF quantiles
%   (default=0.05:0.1:0.95).
%
%   [...,FMODEL] = ANDMODEL(...) returns the AND model based on the joint
%   probability of X and Y (Townsend & Ashby, 1983). By default, the model
%   assumes stochastic independence between RTs on different sensory
%   channels, but this assumption can be specified using the DEP argument
%   (see below). For valid estimates of FMODEL, the stimuli used to
%   generate X, Y and XY should be presented in random order to meet the
%   assumption of context invariance.
%
%   [...,T] = ANDMODEL(...) returns the time intervals used to compute the
%   CDFs for the vertical test.
%
%   [...] = ANDMODEL(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,XY]),max([X,Y,XY])])
%   'dep'       a scalar specifying the model's assumption of stochastic
%               dependence between sensory channels: pass in 0 to assume
%               independence (AND model; default), -1 to assume perfect
%               negative dependence (Colonius-Vorberg lower bound) and 1 to
%               assume perfect positive dependence (Colonius-Vorberg upper
%               bound)
%   'test'      a string specifying how to test the model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%
%   See also ANDMODEL3, ANDGAIN, ANDBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Townsend JT, Ashby FG (1983) Stochastic modeling of elementary
%           psychological processes. Cambridge University Press.
%       [3] Colonius H, Vorberg D (1994) Distribution Inequalities for
%           Parallel Models with Unlimited Capacity. J Math Psychol
%           38:35-58.
%       [4] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 3-May-2019

% Decode input variable arguments
[lim,dep,test] = decode_varargin(varargin);

% Set default values
if ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
elseif nargin < 4 || isempty(p)
    p = 0.05:0.1:0.95;
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
if strcmpi(test,'ver')
    Fx = rt2cdf(x,p,lim);
    Fy = rt2cdf(y,p,lim);
    [Fxy,t] = rt2cdf(xy,p,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fxy = rt2cfp(xy,lim(2));
end

% Compute model
if nargout > 3
    if dep == 0 % AND model
        Fmodel = Fx.*Fy;
    elseif dep == -1 % Colonius-Vorberg lower bound
        Fmodel = max(Fx+Fy-1,zeros(size(Fxy)));
    elseif dep == 1 % Colonius-Vorberg upper bound
        Fmodel = min(Fx,Fy);
    end
    if strcmpi(test,'hor')
        Fmodel = cfp2q(Fmodel,p);
    end
end

% Compute quantiles for horizontal test
if strcmpi(test,'hor')
    Fx = cfp2q(Fx,p);
    Fy = cfp2q(Fy,p);
    Fxy = cfp2q(Fxy,p);
end

% Time intervals for horizontal test not required
if nargout > 4 && strcmpi(test,'hor')
    error('Time intervals T not required for horizontal test.')
end

function [lim,dep,test] = decode_varargin(varargin)
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