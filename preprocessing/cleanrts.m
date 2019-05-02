function [y,idx] = cleanrts(x,cutoffs,varargin)
%cleanrts Clean reaction time data using outlier correction procedures.
%   Y = CLEANRTS(X,CUTOFFS) returns the reaction time data in vector X
%   after the removal of outliers, defined as any values outside of the
%   lower and upper RT cutoffs specified in the 2-element vector CUTOFFS.
%   To not use absolute cutoffs, leave CUTOFFS empty. This function treats
%   NaNs as missing values, and ignores them.
%
%   [...,IDX] = CLEANRTS(...) returns the indices of any element in X
%   identified as an outlier by the chosen procedure.
%
%   [...] = CLEANRTS(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'transform' a string specifying the transformation to apply for the
%               statistical outlier correction method
%                   'none'      no transformation (default)
%                   'inv'       inverse transformation
%                   'log'       log transformation
%   'method'    a string specifying the statistical method to use for
%               outlier correction
%                   'none'      no outlier correction method (default)
%                   'mad'       3 scaled mean absolute deviations (MADs)
%                               above or below the mean
%                   'median'    3 scaled median absolute deviations above
%                               or below the median
%                   'mean'      3 standard deviations above or below the
%                               mean
%                   'iqr'       1.5 interquartile ranges above the upper
%                               quartile or below the lower quartile
%   'factor'    a scalar specifying the factor used to scale the outlier
%               detection threshold specified by METHOD (default=3 or 1.5)
%   'per'       a 2-element vector specifying the lower and upper
%               percentiles of RTs to keep (default=[0,100])
%   'remove'    a scalar specifying whether or not to remove data
%               identified as being outliers: pass in 1 to remove outliers
%               (default) and 0 to keep outliers
%
%   See also RT2CDF, RT2CFP, CFP2PER, GETAUC.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Ratcliff R (1993) Methods for Dealing With Reaction Time
%           Outliers. Psychol Bull 114(3):510-532.
%       [3] Miller J (1991) Reaction Time Analysis with Outlier Exclusion:
%           Bias Varies with Sample Size. Q J Exp Psychol 43A(4):907-912.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 01-May-2019

% Decode input variable arguments
[transform,method,factor,per,remove] = decode_varargin(varargin);

% Transpose if row vector
if isrow(x), x = x'; end

% Identify outliers using cutoffs
if ~isempty(cutoffs)
    idx1 = x<cutoffs(1) | x>cutoffs(2);
else
    idx1 = x<min(x) | x>max(x);
end

% Transform data if specified
if strcmpi(transform,'inv')
    x = 1./x;
elseif strcmpi(transform,'log')
    x = log(x);
end

% Identify outliers using statistical method
if strcmpi(method,'mad')
    mx = nanmean(x);
    madx = nanmean(abs(x-mx));
    idx2 = x<mx-factor*madx | x>mx+factor*madx;
elseif strcmpi(method,'median')
    medx = nanmedian(x);
    madx = nanmedian(abs(x-medx));
    idx2 = x<medx-factor*madx | x>medx+factor*madx;
elseif strcmpi(method,'mean')
    nx = sum(~isnan(x));
    mx = nanmean(x);
    sdx = sqrt(nansum((x-mx).^2)/(nx-1));
    idx2 = x<mx-factor*sdx | x>mx+factor*sdx;
elseif strcmpi(method,'iqr')
    lwr = prctile(x,25);
    upr = prctile(x,75);
    iqr = upr-lwr;
    idx2 = x<lwr-factor*iqr | x>upr+factor*iqr;
else
    idx2 = x<min(x) | x>max(x);
end

% Untransform data
if strcmpi(transform,'inv')
    x = 1./x;
elseif strcmpi(transform,'log')
    x = exp(x);
end

% Identify RTs outside percentile range
per = prctile(x,per);
idx3 = x<per(1) & x>per(2);

% Get indices of outliers
idx = idx1 | idx2 | idx3;

% Remove outliers
if remove
    y = x(idx==0);
else
    y = x;
end

function [transform,method,factor,per,remove] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'transform')) && ~isempty(varargin{find(strcmpi(varargin,'transform'))+1})
    transform = varargin{find(strcmpi(varargin,'transform'))+1};
    if ~any(strcmpi(transform,{'inv','log','none'}))
        error('Invalid value for argument TRANSFORM. Valid values are: ''inv'', ''log'', ''none''.')
    end
else
    transform = 'none'; % default: no transformation
end
if any(strcmpi(varargin,'method')) && ~isempty(varargin{find(strcmpi(varargin,'method'))+1})
    method = varargin{find(strcmpi(varargin,'method'))+1};
    if ~any(strcmpi(method,{'mad','median','mean','iqr'}))
        error('Invalid value for argument METHOD. Valid values are: ''mad'', ''median'', ''mean'', ''iqr''.')
    end
else
    method = 'none'; % default: no outlier correction
end
if any(strcmpi(varargin,'factor')) && ~isempty(varargin{find(strcmpi(varargin,'factor'))+1})
    factor = varargin{find(strcmpi(varargin,'factor'))+1};
    if ~isnumeric(factor) || ~isscalar(factor) || isnan(factor) || isinf(factor) || factor<=0
        error('FACTOR must be a scalar of positive value.')
    end
else
    if any(strcmpi(method,{'mad','median','mean'}))
        factor = 3; % default: 3
    elseif strcmpi(method,{'iqr'})
        factor = 1.5; % default: 1.5
    else
        factor = []; % default: []
    end
end
if any(strcmpi(varargin,'per')) && ~isempty(varargin{find(strcmpi(varargin,'per'))+1})
    per = varargin{find(strcmpi(varargin,'per'))+1};
    if ~isnumeric(per) || isscalar(per) || any(isnan(per)) || any(isinf(per)) || any(per<0) || any(per>100) || per(1)>=per(2)
        error('PER must be a 2-element vector with values between 0 and 100.')
    end
else
    per = [0,100]; % default: keep all RTs
end
if any(strcmpi(varargin,'remove')) && ~isempty(varargin{find(strcmpi(varargin,'remove'))+1})
    remove = varargin{find(strcmpi(varargin,'remove'))+1};
    if remove~=0 && remove~=1
        error('REMOVE must be a scalar with a value of 0 or 1.')
    end
else
    remove = 1; % default: remove outliers
end