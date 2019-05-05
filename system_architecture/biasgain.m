function [gain,Fdiff,t] = biasgain(Xx,Xy,Xxy,Yx,Yy,Yxy,XY,p,varargin)
%biasgain Multisensory gain with a bias towards X or Y.
%   GAIN = BIASGAIN(XX,XY,XXY,YX,YY,YXY,XY) returns the multisensory gain
%   for a bisensory OR or AND task, quantified by the area between the CDFs
%   of the bisensory RT distribution XY, and the bias model based on the
%   mean probability of XX, XY and XXY (Crosse et al., 2019). By default,
%   the model is biased towards X, but this bias can be specified using the
%   BIAS argument (see below). X, Y and XY can have different lengths. This
%   function treats NaNs as missing values, and ignores them.
%
%   [...] = BIASGAIN(...,P) uses the intervals P to generate CDFs. P is a
%   vector of decimal values between 0 and 1 inclusive. For horizontal
%   tests, P is the probabilities used to compute the CDF quantiles
%   (default=0.05:0.1:0.95).
%
%   [...,FDIFF] = BIASGAIN(...) returns the difference at each quantile to
%   test for violations of the model.
%
%   [...,T] = BIASGAIN(...) returns the time intervals used to compute the
%   CDFs for the vertical test.
%
%   [...] = BIASGAIN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly to other conditions
%               (default=[min([X,Y,XY]),max([X,Y,XY])])
%   'bias'      a string specifying how to bias the model: 1) bias towards
%               a specific (dominant) modality (X or Y), or 2) bias towards
%               the previous modality (n-1), except when the previous
%               modality is XY (bias towards X or Y)
%                   'X'         X bias (default)
%                   'Y'         Y bias
%                   'n-1X'      n-1 bias (X bias when n-1=XY)
%                   'n-1Y'      n-1 bias (Y bias when n-1=XY)
%   'test'      a string specifying how to test the competition model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%
%   See also TRIALHISTORY, BIASMODEL, BIASBENEFIT, TPERMTEST, EFFECTSIZE.
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
[lim,bias,test] = decode_varargin(varargin);

% Set default values
if ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
elseif nargin < 8 || isempty(p)
    p = 0.05:0.1:0.95;
end

% Transpose row vectors
if isrow(Xx), Xx = Xx'; end
if isrow(Xy), Xy = Xy'; end
if isrow(Xxy), Xxy = Xxy'; end
if isrow(Yx), Yx = Yx'; end
if isrow(Yy), Yy = Yy'; end
if isrow(Yxy), Yxy = Yxy'; end
if isrow(XY), XY = XY'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([Xx;Xy;Xxy;Yx;Yy;Yxy;XY]),max([Xx;Xy;Xxy;Yx;Yy;Yxy;XY])];
end

% Compute CDFs
if strcmpi(test,'ver')
    FXx = rt2cdf(Xx,p,lim);
    FXy = rt2cdf(Xy,p,lim);
    FXxy = rt2cdf(Xxy,p,lim);
    FYx = rt2cdf(Yx,p,lim);
    FYy = rt2cdf(Yy,p,lim);
    FYxy = rt2cdf(Yxy,p,lim);
    [FXY,t] = rt2cdf(XY,p,lim);
elseif strcmpi(test,'hor')
    FXx = rt2cfp(Xx,lim(2));
    FXy = rt2cfp(Xy,lim(2));
    FXxy = rt2cfp(Xxy,lim(2));
    FYx = rt2cfp(Yx,lim(2));
    FYy = rt2cfp(Yy,lim(2));
    FYxy = rt2cfp(Yxy,lim(2));
    FXY = rt2cfp(XY,lim(2));
end

% Compute model
if strcmpi(bias,'X') % X bias
    Fmodel = (FXx+FXy+FXxy)/3;
elseif strcmpi(bias,'Y') % Y bias
    Fmodel = (FYx+FYy+FYxy)/3;
elseif strcmpi(bias,'n-1X') % n-1 bias (X bias when n-1=XY)
    Fmodel = (FXx+FYy+FXxy)/3;
elseif strcmpi(bias,'n-1Y') % n-1 bias (Y bias when n-1=XY)
    Fmodel = (FXx+FYy+FYxy)/3;
end

% Compute quantiles for horizontal test
if strcmpi(test,'hor')
    FXY = cfp2q(Fmodel,p);
    Fmodel = cfp2q(Fmodel,p);
end

% Compute difference
if strcmpi(test,'ver')
    Fdiff = FXY-Fmodel;
elseif strcmpi(test,'hor')
    Fdiff = Fmodel-FXY;
end

% Compute multisensory gain
gain = getauc(p,Fdiff,area);

% Time intervals for horizontal test not required
if nargout > 2 && strcmpi(test,'hor')
    error('Time intervals T not required for horizontal test.')
end

function [lim,bias,test] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'bias')) && ~isempty(varargin{find(strcmpi(varargin,'bias'))+1})
    bias = varargin{find(strcmpi(varargin,'bias'))+1};
    if ~any(strcmpi(bias,{'X','Y','n-1X','n-1Y'}))
        error('Invalid value for argument BIAS. Valid values are: ''X'', ''Y'', ''n-1X'', ''n-1Y''.')
    end
else
    bias = 'X'; % default: X bias
end
if any(strcmpi(varargin,'test')) && ~isempty(varargin{find(strcmpi(varargin,'test'))+1})
    test = varargin{find(strcmpi(varargin,'test'))+1};
    if ~any(strcmpi(test,{'ver','hor'}))
        error('Invalid value for argument TEST. Valid values are: ''ver'', ''hor''.')
    end
else
    test = 'ver'; % default: vertical test
end