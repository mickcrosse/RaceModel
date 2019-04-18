function [MIC,SIC,q,lim] = sft(xhh,xll,xhl,xlh,varargin)
%sft Systems factorial technology.
%   MIC = SFT(XHH,XLL,XHL,XLH) returns the mean interaction contrast (MIC)
%   based on the mean values of the bisensory conditions XHH, XLL, XHL and
%   XLH, where H indicates high salience and L indicates low salience on
%   the respective channels. XHH, XLL, XHL and XLH are not required to have
%   an equal number of observations. This function treats NaNs as missing
%   values, and ignores them.
%
%   [...,SIC] = SFT(...) returns survivor interaction contrast (SIC) based
%   on the survivor functions of the bisensory conditions XHH, XLL, XHL and
%   XLH.
%
%   The architecture and stopping rule can be inferred by the combined
%   outcome of MIC and SIC using the following reference table:
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
%   'outlier'   a 2-element vector specifying the lower and upper RT
%               cutoffs for outlier correction (default=no correction).
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
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
%           and coactive theories. J Math Psychol 39(4):321�359.
%       [3] Little DR, Altieri N, Fific M, Yang CT (2017) Systems factorial
%           technology: A theory driven methodology for the identification
%           of perceptual and cognitive mechanisms. Academic Press.
%
%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 18-Apr-2019

% Decode input variable arguments
[p,outlier,per,lim] = decode_varargin(varargin);

% Outlier correction procedure
if ~isempty(outlier)
    xhh(xhh<outlier(1)|xhh>outlier(2)) = [];
    xhl(xhl<outlier(1)|xhl>outlier(2)) = [];
    xlh(xlh<outlier(1)|xlh>outlier(2)) = [];
    xll(xll<outlier(1)|xll>outlier(2)) = [];
end

% Get RT range for each condition
lims = zeros(4,2);
lims(1,:) = prctile(xhh,per);
lims(2,:) = prctile(xhl,per);
lims(3,:) = prctile(xlh,per);
lims(4,:) = prctile(xll,per);

% Limit RTs to specified range
xhh = xhh(xhh>=lims(1,1) & xhh<=lims(1,2));
xhl = xhl(xhl>=lims(2,1) & xhl<=lims(2,2));
xlh = xlh(xlh>=lims(3,1) & xlh<=lims(3,2));
xll = xll(xll>=lims(4,1) & xll<=lims(4,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
Fxhh = rt2cdf(xhh,p,lim);
Fxhl = rt2cdf(xhl,p,lim);
Flh = rt2cdf(xlh,p,lim);
[Fxll,q] = rt2cdf(xll,p,lim);

% Compute mean interaction contrast
MIC = nanmean(xhh)-nanmean(xhl)-nanmean(xlh)+nanmean(xll);

% Compute survivor interaction contrast
SIC = Fxhl+Flh-Fxhh-Fxll;

function [p,outlier,per,lim] = decode_varargin(varargin)
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