function [gain,Fdiff,q,lim] = orgain(x,y,xy,varargin)
%orgain Multisensory gain for a bisensory OR task.
%   GAIN = ORGAIN(X,Y,XY) returns the multisensory gain for a bisensory OR
%   task, quantified by the area between the CDFs of the bisensory RT
%   distribution XY, and the OR (race) model based on the probability
%   summation of X and Y (Colonius & Diederich, 2006). X, Y and XY can have
%   different lengths. This function treats NaNs as missing values, and
%   ignores them.
%
%   [...,FDIFF] = ORGAIN(...) returns the difference at each quantile to
%   test for violations of the model.
%
%   [...,Q] = ORGAIN(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test or the probabilities used to compute the
%   percentiles for the horizontal test.
%
%   [...,LIM] = ORGAIN(...) returns the lower and upper RT limits used to
%   compute the CDFs. These values can be used to set the CDF limits of
%   subsequent tests that are to be compared with this one.
%
%   [...] = ORGAIN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%   'dep'       a scalar specifying the model's assumption of stochastic
%               dependence between sensory channels: pass in 0 to assume
%               independence (OR model; default), -1 to assume perfect
%               negative dependence (Miller's bound) and 1 to assume
%               perfect positive dependence (Grice's bound)
%   'test'      a string specifying how to test the model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%
%   See also ORGAIN3, ORMODEL, ORBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Colonius H, Diederich A (2006) The Race Model Inequality:
%           Interpreting a Geometric Measure of the Amount of Violation.
%           Psychol Rev 113(1):148�154.
%       [3] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 01-May-2019

% Decode input variable arguments
[p,lim,dep,test,area] = decode_varargin(varargin);

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
    [Fxy,q] = rt2cdf(xy,p,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fxy = rt2cfp(xy,lim(2));
end

% Compute model
if dep == 0 % OR model
    Fmodel = Fx+Fy-Fx.*Fy;
elseif dep == -1 % Miller's bound
    Fmodel = min(Fx+Fy,ones(size(Fxy)));
elseif dep == 1 % Grice's bound
    Fmodel = max(Fx,Fy);
end

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    Fxy = cfp2per(Fxy,p);
    Fmodel = cfp2per(Fmodel,p);
end

% Compute difference
if strcmpi(test,'ver')
    Fdiff = Fxy-Fmodel;
elseif strcmpi(test,'hor')
    Fdiff = Fmodel-Fxy;
end

% Compute multisensory gain
gain = getauc(p,Fdiff,area);

% Get y-values for horizontal test
if nargout > 2 &&  strcmpi(test,'hor')
    q = p;
end

function [p,lim,dep,test,area] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'dep')) && ~isempty(varargin{find(strcmpi(varargin,'dep'))+1})
    dep = varargin{find(strcmpi(varargin,'dep'))+1};
    if dep~=0 && dep~=-1 && dep~=1
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