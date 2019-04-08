function [MSE,Fdiff,q] = switchcost(sw,re,varargin)
%switchcost Modality switch effect of RTs for mixed stimulus presenation.
%   MSE = SWITCHCOST(SW,RE) returns the modality switch effect (MSE) of
%   response times for mixed stimulus presenation. MSE is quantified as the
%   area between the cumulative distribution functions (CDFs) of the RT
%   distributions of the switch trials SW, and repeat trials RE (Crosse et
%   al., 2019). SW and RE are not required to have an equal number of
%   observations. This function treats NaNs as missing values, and ignores
%   them.
%
%   FDIFF = SWITCHCOST(...) returns the difference between the CDFs of the
%   switch and repeat trials at every quantile.
%
%   [...,Q] = SWITCHCOST(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test.
%
%   [...] = SWITCHCOST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
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
%               unspecified unless comparing directly to other conditions
%               (default=[min([SW,RE]),max([SW,RE])])
%   'test'      a string specifying how to test the MSE
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%
%   See also RACEMODEL, RSEGAIN, RSEBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) Developmental Recovery of
%           Impaired Multisensory Processing in Autism and the Cost of
%           Switching Sensory Modality. bioRxiv 10.1101/565333.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 4-Apr-2019

% Decode input variable arguments
[p,outlier,per,lim,test,area] = decode_varargin(varargin);

% Outlier correction procedure
if ~isempty(outlier)
    sw(sw<outlier(1)|sw>outlier(2)) = [];
    re(re<outlier(1)|re>outlier(2)) = [];
end

% Get RT range for each condition
lims = zeros(3,2);
lims(1,:) = prctile(sw,per);
lims(2,:) = prctile(re,per);

% Limit RTs to specified range
sw = sw(sw>=lims(1,1) & sw<=lims(1,2));
re = re(re>=lims(2,1) & re<=lims(2,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
if strcmpi(test,'ver')
    Fsw = rt2cdf(sw,p,lim);
    [Fre,q] = rt2cdf(re,p,lim);
elseif strcmpi(test,'hor')
    Fsw = rt2cfp(sw,lim(2));
    Fre = rt2cfp(re,lim(2));
end

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    Fsw = cfp2per(Fsw,p);
    Fre = cfp2per(Fre,p);
end

% Compute difference
if strcmpi(test,'ver')
    Fdiff = Fre-Fsw;
elseif strcmpi(test,'hor')
    Fdiff = Fsw-Fre;
end

% Compute MSE
MSE = getauc(p,Fdiff,area);

function [p,outlier,per,lim,test,area] = decode_varargin(varargin)
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