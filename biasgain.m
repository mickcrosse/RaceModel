function [gain,Fdiff,q,lim] = biasgain(Xx,Xy,Xxy,Yx,Yy,Yxy,XY,varargin)
%biasgain Multisensory gain with a bias towards X or Y.
%   GAIN = BIASGAIN(XX,XY,XXY,YX,YY,YXY,XY) returns the multisensory gain
%   for a bisensory OR or AND task, quantified by the area between the CDFs
%   of the bisensory RT distribution XY, and the bias model based on the
%   mean probability of XX, XY and XXY (Crosse et al., 2019). By default,
%   the model is biased towards X, but this bias can be specified using the
%   BIAS argument (see below). X, Y and XY are not required to have an
%   equal number of observations. This function treats NaNs as missing
%   values, and ignores them.
%
%   [...,FDIFF] = BIASGAIN(...) returns the difference at each quantile to
%   test for violations of the model.
%
%   [...,Q] = BIASGAIN(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test and the probabilities used to compute the
%   percentiles for the horizontal test.
%
%   [...,LIM] = BIASGAIN(...) returns the lower and upper RT limits used to
%   compute the CDFs. These values can be used to set the CDF limits of
%   subsequent tests that are to be compared with this one.
%
%   [...] = BIASGAIN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%   Apr 2017; Last Revision: 14-Apr-2019

% Decode input variable arguments
[p,outlier,per,lim,bias,test] = decode_varargin(varargin);

% Outlier correction procedure
if ~isempty(outlier)
    Xx(Xx<outlier(1)|Xx>outlier(2)) = [];
    Xy(Xy<outlier(1)|Xy>outlier(2)) = [];
    Xxy(Xxy<outlier(1)|Xxy>outlier(2)) = [];
    Yx(Yx<outlier(1)|Yx>outlier(2)) = [];
    Yy(Yy<outlier(1)|Yy>outlier(2)) = [];
    Yxy(Yxy<outlier(1)|Yxy>outlier(2)) = [];
    XY(XY<outlier(1)|XY>outlier(2)) = [];
end

% Get RT range for each condition
lims = zeros(7,2);
lims(1,:) = prctile(Xx,per);
lims(2,:) = prctile(Xy,per);
lims(3,:) = prctile(Xxy,per);
lims(4,:) = prctile(Yx,per);
lims(5,:) = prctile(Yy,per);
lims(6,:) = prctile(Yxy,per);
lims(7,:) = prctile(XY,per);

% Limit RTs to specified range
Xx = Xx(Xx>=lims(1,1) & Xx<=lims(1,2));
Xy = Xy(Xy>=lims(2,1) & Xy<=lims(2,2));
Xxy = Xxy(Xxy>=lims(3,1) & Xxy<=lims(3,2));
Yx = Yx(Yx>=lims(4,1) & Yx<=lims(4,2));
Yy = Yy(Yy>=lims(5,1) & Yy<=lims(5,2));
Yxy = Yxy(Yxy>=lims(6,1) & Yxy<=lims(6,2));
XY = Xy(Xy>=lims(7,1) & XY<=lims(7,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
if strcmpi(test,'ver')
    FXx = rt2cdf(Xx,p,lim);
    FXy = rt2cdf(Xy,p,lim);
    FXxy = rt2cdf(Xxy,p,lim);
    FYx = rt2cdf(Yx,p,lim);
    FYy = rt2cdf(Yy,p,lim);
    FYxy = rt2cdf(Yxy,p,lim);
    [FXY,q] = rt2cdf(XY,p,lim);
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

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    FXY = cfp2per(Fmodel,p);
    Fmodel = cfp2per(Fmodel,p);
end

% Compute difference
if strcmpi(test,'ver')
    Fdiff = FXY-Fmodel;
elseif strcmpi(test,'hor')
    Fdiff = Fmodel-FXY;
end

% Compute multisensory gain
gain = getauc(p,Fdiff,area);

% Get y-values for horizontal test
if nargout > 2 &&  strcmpi(test,'hor')
    q = p;
end

function [p,outlier,per,lim,bias,test] = decode_varargin(varargin)
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