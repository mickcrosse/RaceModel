function [MIC,SIC,q,lim] = sft(xhh,xll,xhl,xlh,varargin)
%sft Systems factorial technology.
%   MIC = SFT(XHH,XLL,XHL,XLH) returns the mean interaction contrast (MIC)
%   based on the mean values of the bisensory conditions XHH, XLL, XHL and
%   XLH, where H indicates high salience and L indicates low salience on
%   the respective channels. XHH, XLL, XHL and XLH can have different
%   lengths. This function treats NaNs as missing values, and ignores them.
%
%   [...,SIC] = SFT(...) returns survivor interaction contrast (SIC) based
%   on the survivor functions of the bisensory conditions XHH, XLL, XHL and
%   XLH.
%
%   The system's architecture and stopping rule can be inferred from the
%   following reference table using the combined outcome of MIC and SIC.
%
%   Model           MIC     SIC
%   Serial, OR      0       0
%   Serial, AND     0       -to+
%   Parallel, OR	>0      >0
%   Parallel, AND	<0      <0
%   Coactive        >0      -to+
%
%   [...,Q] = SFT(...) returns the RT quantiles used to compute the CDFs.
%
%   [...,LIM] = SFT(...) returns the lower and upper RT limits used to
%   compute the CDFs. These values can be used to set the CDF limits of
%   subsequent tests that are to be compared with this one.
%
%   [...] = SFT(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values. Valid parameters are the following:
%
%   Parameter   Value
%   'p'         a vector specifying the probabilities for computing the
%               quantiles of a vertical test or the percentiles of a
%               horizontal test (default=0.05:0.1:0.95)
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'verbose'   a scalar specifying whether to display the reference table
%               for inferring the system's architecture and stopping rule:
%               pass in 1 to display reference table (default) and 0 to not
%
%   See also ORCAPACITY, ANDCAPACITY, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Townsend JT, Nozawa G (1995) Spatio-temporal properties of
%           elementary perception: An investigation of parallel, serial,
%           and coactive theories. J Math Psychol 39(4):321–359.
%       [3] Little DR, Altieri N, Fific M, Yang CT (2017) Systems factorial
%           technology: A theory driven methodology for the identification
%           of perceptual and cognitive mechanisms. Academic Press.
%
%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 01-May-2019

% Decode input variable arguments
[p,lim,verbose] = decode_varargin(varargin);

% Generate reference table
Model = {'Serial, OR';'Serial, AND';'Parallel, OR';'Parallel, AND';'Coactive'};
MIC = {'0';'0';'>0';'<0';'>0'};
SIC = {'0';'-to+';'>0';'<0';'-to+'};
reftable = table(Model,MIC,SIC);

% Transpose row vectors
if isrow(xhh), xhh = xhh'; end
if isrow(xll), xll = xll'; end
if isrow(xhl), xhl = xhl'; end
if isrow(xlh), xlh = xlh'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([xhh;xll;xhl;xlh]),max([xhh;xll;xhl;xlh])];
end

% Compute CDFs
Fxhh = rt2cdf(xhh,p,lim);
Fxhl = rt2cdf(xhl,p,lim);
Flh = rt2cdf(xlh,p,lim);
[Fxll,q] = rt2cdf(xll,p,lim);

% Compute mean interaction contrast
MIC = nanmean(xhh)+nanmean(xll)-nanmean(xhl)-nanmean(xlh);

% Compute survivor interaction contrast
SIC = Fxhl+Flh-Fxhh-Fxll;

% Display reference table
if verbose
    disp(reftable)
end

function [p,lim,verbose] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'verbose')) && ~isempty(varargin{find(strcmpi(varargin,'verbose'))+1})
    verbose = varargin{find(strcmpi(varargin,'verbose'))+1};
    if verbose~=0 && verbose~=1
        error('VERBOSE must be a scalar with a value of 0 or 1.')
    end
else
    verbose = 0; % default: display table
end