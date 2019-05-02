function [Bemp,Bpred,Femp,Fpred,q,lim] = biasbenefit(Xx,Xy,Xxy,Yx,Yy,Yxy,XY,varargin)
%biasbenefit Multisensory benefit with a bias towards X or Y.
%   BEMP = BIASBENEFIT(XX,XY,XXY,YX,YY,YXY,XY) returns the empirical
%   multisensory benefit for a bisensory OR task, quantified by the area
%   between the CDFs of the faster of the unisensory RT distributions X and
%   Y, and the bisensory RT distribution XY (Otto et al., 2013). By
%   default, the bound used to compute benefits is for an OR task design,
%   but the task design can be specified using the TASK argument (see
%   below). X, Y and XY can have different lengths. This function treats
%   NaNs as missing values, and ignores them.
%
%   [...,BPRED] = BIASBENEFIT(...) returns the predicted multisensory
%   benefit for a bisensory OR task, quantified by the area between the
%   CDFs of the faster of the unisensory RT distributions X and Y, and the
%   bias model based on the mean probability of XX, XY and XXY (Crosse et
%   al., 2019a,b). By default, the model is biased towards X, but this bias
%   can be specified using the BIAS argument (see below).
%
%   [...,FEMP,FPRED] = BIASBENEFIT(...) returns the difference at each
%   quantile for empirical and predicted benefits, respectively.
%
%   [...,Q] = BIASBENEFIT(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test and the probabilities used to compute the
%   percentiles for the horizontal test.
%
%   [...,LIM] = BIASBENEFIT(...) returns the lower and upper RT limits used
%   to compute the CDFs. These values can be used to set the CDF limits of
%   subsequent tests that are to be compared with this one.
%
%   [...] = BIASBENEFIT(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'p'         a vector specifying the probabilities for computing the
%               quantiles of a vertical test or the percentiles of a
%               horizontal test (default=0.05:0.1:0.95)
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
%   'task'      a string specifying the task design
%                   'OR'        OR task (default)
%                   'AND'       AND task
%
%   See also TRIALHISTORY, BIASMODEL, BIASGAIN, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019a) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Crosse MJ, Foxe JJ, Molholm S (2019b) Developmental Recovery of
%           Impaired Multisensory Processing in Autism and the Cost of
%           Switching Sensory Modality. bioRxiv 10.1101/565333.
%       [3] Otto TU, Dassy B, Mamassian P (2013) Principles of multisensory
%           behavior. J Neurosci 33(17):7463-7474.
%       [4] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 01-May-2019

% Decode input variable arguments
[p,lim,bias,test,task] = decode_varargin(varargin);

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

% Compute bound
Fx = mean([FXx,FXy,FXxy],2);
Fy = mean([FYx,FYy,FYxy],2);
if strcmpi(task,'OR') % Grice's bound
    Fmax = max(Fx,Fy);
elseif strcmpi(task,'AND') % Colonius-Vorberg upper bound
    Fmax = min(Fx,Fy);
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
    FXY = cfp2per(FXY,p);
    Fmax = cfp2per(Fmax,p);
    Fmodel = cfp2per(Fmodel,p);
end

% Compute difference
if strcmpi(test,'ver')
    Femp = FXY-Fmax;
elseif strcmpi(test,'hor')
    Femp = Fmax-FXY;
end

% Compute empirical benefit
Bemp = getauc(p,Femp,area);

if nargout > 1
    
    % Compute difference
    if strcmpi(test,'ver')
        Fpred = Fmodel-Fmax;
    elseif strcmpi(test,'hor')
        Fpred = Fmax-Fmodel;
    end
    
    % Compute predicted benefit
    Bpred = getauc(p,Fpred,area);
    
end

% Get y-values for horizontal test
if nargout > 4 &&  strcmpi(test,'hor')
    q = p;
end

function [p,lim,bias,test,task] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'task')) && ~isempty(varargin{find(strcmpi(varargin,'task'))+1})
    task = varargin{find(strcmpi(varargin,'task'))+1};
    if ~any(strcmpi(task,{'OR','AND'}))
        error('Invalid value for argument TASK. Valid values are: ''OR'', ''AND''.')
    end
else
    task = 'ver'; % default: OR task
end