function bpred = competemodel(Xx,Xy,Xxy,Yx,Yy,Yxy,varargin)
%competemodel Multisensory benefit predicted by competition models.
%   BPRED = COMPETEMODEL(XX,XY,XXY,YX,YY,YXY) returns the predicted benefit
%   of multisensory competition, quantified by the area between the CDFs of
%   the most effective of the unisensory RT distributions X and Y, and the
%   competition model biased towards X. Unisensory RT distributions X and Y
%   should be separated by previous modality: XX, XY, XXY and YX, YY, YXY,
%   respectively. Competition can be modelled as a bias towards a specific
%   (i.e., dominant) modality (X or Y) or a bias towards the previous
%   modality (n-1), except when the previous modality is XY (biased towards
%   either X or Y). For a mathematical description, see Crosse et al.
%   (2019). XX, XY, XXY, YX, YY and YXY are not required to have an equal
%   number of observations. This function treats NaNs as missing values,
%   and ignores them.
%
%   [...] = COMPETEMODEL(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%               (default=[min([X,Y]),max([X,Y])])
%   'dep'       a scalar specifying the model's assumption of statistical
%               dependence between sensory channels: pass in 0 to assume
%               independence (Raab, 1962; default) and -1 to assume a
%               perfect negative dependence (Miller, 1982)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%   'model'     a string specifying whether competition is modelled as 1)
%               a bias towards a specific (dominant) modality (X or Y), or
%               2) a bias towards the previous modality (n-1), except when
%               the previous modality is XY (biased towards either X or Y)
%                   '1X'        X bias for all trials (default)
%                   '1Y'        Y bias for all trials
%                   '2X'        n-1 bias and X bias when previous is XY
%                   '2Y'        n-1 bias and Y bias when previous is XY
%
%   See also COMPETEMODEL3, RSEBENEFIT, RSEGAIN, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) Developmental Recovery of
%           Impaired Multisensory Processing in Autism and the Cost of
%           Switching Sensory Modality. bioRxiv 10.1101/565333.
%       [2] Otto TU, Dassy B, Mamassian P (2013) Principles of multisensory
%           behavior. J Neurosci 33(17):7463-7474.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 4-Apr-2019

% Decode input variable arguments
[p,outlier,per,lim,test,area] = decode_varargin(varargin);

% Outlier correction procedure
if ~isempty(outlier)
    Xx(Xx<outlier(1)|Xx>outlier(2)) = [];
    Xy(Xy<outlier(1)|Xy>outlier(2)) = [];
    Xxy(Xxy<outlier(1)|Xxy>outlier(2)) = [];
    Yx(Yx<outlier(1)|Yx>outlier(2)) = [];
    Yy(Yy<outlier(1)|Yy>outlier(2)) = [];
    Yxy(Yxy<outlier(1)|Yxy>outlier(2)) = [];
end

% Get RT range for each condition
lims = zeros(6,2);
lims(1,:) = prctile(Xx,per);
lims(2,:) = prctile(Xy,per);
lims(3,:) = prctile(Xxy,per);
lims(4,:) = prctile(Yx,per);
lims(5,:) = prctile(Yy,per);
lims(6,:) = prctile(Yxy,per);

% Limit RTs to specified range
Xx = Xx(Xx>=lims(1,1) & Xx<=lims(1,2));
Xy = Xy(Xy>=lims(2,1) & Xy<=lims(2,2));
Xxy = Xxy(Xxy>=lims(3,1) & Xxy<=lims(3,2));
Yx = Yx(Yx>=lims(4,1) & Yx<=lims(4,2));
Yy = Yy(Yy>=lims(5,1) & Yy<=lims(5,2));
Yxy = Yxy(Yxy>=lims(6,1) & Yxy<=lims(6,2));

% Get pooled distributions
if size(Xx,2) == 1
    x = [Xx;Xy;Xxy];
    y = [Yx;Yy;Yxy];
elseif size(Xx,1) == 1
    x = [Xx,Xy,Xxy];
    y = [Yx,Yy,Yxy];
end

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
if strcmpi(test,'ver')
    Fx = rt2cdf(x,p,lim);
    Fy = rt2cdf(y,p,lim);
    Fxx = rt2cdf(Xx,p,lim);
    Fxy = rt2cdf(Xy,p,lim);
    Fxxy = rt2cdf(Xxy,p,lim);
    Fyx = rt2cdf(Yx,p,lim);
    Fyy = rt2cdf(Yy,p,lim);
    Fyxy = rt2cdf(Yxy,p,lim);
elseif strcmpi(test,'hor')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fxx = rt2cfp(Xx,lim(2));
    Fxy = rt2cfp(Xy,lim(2));
    Fxxy = rt2cfp(Xxy,lim(2));
    Fyx = rt2cfp(Yx,lim(2));
    Fyy = rt2cfp(Yy,lim(2));
    Fyxy = rt2cfp(Yxy,lim(2));
end

% Compute Grice's bound
Fmax = max(Fx,Fy);

% Compute competition model
if strcmpi(bias,'1X')
    Fcomp = (Fxx+Fxy+Fxxy)/3;
elseif strcmpi(bias,'1Y')
    Fcomp = (Fyx+Fyy+Fyxy)/3;
elseif strcmpi(bias,'2X')
    Fcomp = (Fxx+Fyy+Fxxy)/3;
elseif strcmpi(bias,'2Y')
    Fcomp = (Fxx+Fyy+Fyxy)/3;
end

% Compute percentiles for horizontal test
if strcmpi(test,'hor')
    Fmax = cfp2per(Fmax,p);
    Fcomp = cfp2per(Fcomp,p);
end

% Compute difference
if strcmpi(test,'ver')
    Fdiff = Fcomp-Fmax;
elseif strcmpi(test,'hor')
    Fdiff = Fmax-Fcomp;
end

% Compute predicted benefit
bpred = getauc(p,Fdiff,area);

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