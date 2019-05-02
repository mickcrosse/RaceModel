function [cost,Fdiff,q] = switchcost(sw,re,varargin)
%switchcost Switch cost of reaction times to mixed stimuli.
%   COST = SWITCHCOST(SW,RE) returns the switch cost of reaction times to
%   mixed stimuli, quantified as the area between the CDFs of the RT
%   distributions of the switch trials SW, and repeat trials RE (Crosse et
%   al., 2019a,b). SW and RE can have different lengths. This function
%   treats NaNs as missing values, and ignores them.
%
%   FDIFF = SWITCHCOST(...) returns the difference between the CDFs of the
%   switch and repeat trials at every quantile.
%
%   [...,Q] = SWITCHCOST(...) returns the RT quantiles used to compute the
%   CDFs for the vertical test and the probabilities used to compute the
%   percentiles for the horizontal test.
%
%   [...] = SWITCHCOST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%               (default=[min([SW,RE]),max([SW,RE])])
%   'test'      a string specifying how to test the MSE
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test (Ulrich et al., 2007)
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%
%   See also TRIALHISTORY, ORGAIN, BIASMODEL, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019a) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Crosse MJ, Foxe JJ, Molholm S (2019b) Developmental Recovery of
%           Impaired Multisensory Processing in Autism and the Cost of
%           Switching Sensory Modality. bioRxiv 10.1101/565333.
%       [3] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 01-May-2019

% Decode input variable arguments
[p,lim,test,area] = decode_varargin(varargin);

% Transpose row vectors
if isrow(sw), sw = sw'; end
if isrow(re), re = re'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([sw;re]),max([sw;re])];
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
cost = getauc(p,Fdiff,area);

% Get probabilities for horizontal test
if nargout > 2 &&  strcmpi(test,'hor')
    q = p;
end

function [p,lim,test,area] = decode_varargin(varargin)
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