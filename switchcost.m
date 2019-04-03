function [mse,fdiff] = switchcost(sw,re,varargin)
%switchcost Modality switch effect for mixed multisensory stimuli.
%   MSE = SWITCHCOST(SW,RE) returns the modality switch effect (MSE) of RTs
%   to mixed multisensory stimuli, quantified by the area between the
%   cumulative distribution functions (CDFs) of the RT distributions of the
%   swtich trials SW, and repeat trials RE (Crosse et al., 2019). This
%   function does not require SW and RE to have an equal number of
%   observations. This function treats NaNs as missing values, and ignores
%   them.
%
%   FDIFF = SWITCHCOST(...) returns the difference between the CDFs of the
%   switch and repeat trials at every quantile.
%
%   [...] = SWITCHCOST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'q'         a vector specifying the quantiles to be used to compute the
%               CDFs (default=[0.05:0.05:1])
%   'per'       a 2-element vector specifying the lower and upper RT
%               percentiles to be used for each condition (default=[0,100])
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               to be used to compute the CDFs: it is recommended to leave
%               this unspecified unless comparing directly to other
%               conditions (default=[min([sw,re]),max([sw,re])])
%   'test'      a string specifying how to test the race model
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
%   Apr 2017; Last Revision: 3-Apr-2019

% Decode input variable arguments
[q,per,lim,test,area] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(3,2);
lims(1,:) = prctile(sw,per);
lims(2,:) = prctile(re,per);

% Limit RTs to specified range
sw = sw(sw>lims(1,1) & sw<lims(1,2));
re = re(re>lims(2,1) & re<lims(2,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
fsw = rt2cdf(sw,q,lim);
fre = rt2cdf(re,q,lim);

% Compute difference
if strcmpi(test,'ver')
    fdiff = fre-fsw;
elseif strcmpi(test,'hor')
    fdiff = fsw-fre;
end

% Compute MSE
mse = getauc(q,fdiff,area);

function [q,per,lim,dep,test,area] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'q')) && ~isempty(varargin{find(strcmpi(varargin,'q'))+1})
    q = varargin{find(strcmpi(varargin,'q'))+1};
    if ~isnumeric(q) || isscalar(q) || any(isnan(q)) || any(isinf(q)) || any(q<0) || any(q>1) || any(diff(q)<=0)
        error('Q must be a vector with values between 0 and 1.')
    end
else
    q = 0.05:0.05:1; % default: 0.05 to 1 in 0.05 increments
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
    elseif strcmpi(test,'hor')
        error('Horizontal test not yet implemented. Please watch out for updates.')
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