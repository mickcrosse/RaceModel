function [bemp,bpred] = rsebenefit3(x,y,z,xyz,varargin)
%rsebenefit3 Multisensory benefit of a trisensory redundant signals effect.
%   BEMP = RSEBENEFIT3(X,Y,Z,XYZ) returns the empirical benefit of a
%   redundant signals effect (RSE), quantified by the area between the
%   cumulative distribution functions (CDFs) of the most effective of the
%   unisensory RT distributions X, Y and Z, and the trisensory RT
%   distribution XYZ (Otto et al., 2013). This function does not require X,
%   Y, Z and XYZ to have an equal number of observations. This function
%   treats NaNs as missing values, and ignores them.
%
%   To compute the benefit for the bisensory conditions XY, XZ and YZ, use
%   the function RSEBENEFIT on the corresponding unisensory and bisensory
%   RTs. To compare across bisensory and trisensory conditions, use the
%   same lower and upper RT limits (see below).
%
%   [...,BPRED] = RSEBENEFIT3(...) returns the predicted benefit of an RSE,
%   quantified by the area between the CDFs of the most effective of the
%   unisensory RT distributions X, Y and Z, and the trisensory race model
%   based on X, Y and Z (Otto et al., 2013).
%
%   [...] = RSEBENEFIT3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%               conditions (default=[min([x,y,z,xyz]),max([x,y,z,xyz])])
%   'dep'       a scalar specifying the model's assumption of statistical
%               dependence between X and Y: pass in 0 to assume
%               independence (Raab, 1962; default), -1 to assume a perfect
%               negative dependence (Miller, 1982) and 1 to assume a
%               perfect positive dependence (Grice et al., 1986)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%
%   See also RSEBENEFIT, RACEMODEL3, RSEGAIN3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Otto TU, Dassy B, Mamassian P (2013) Principles of multisensory
%           behavior. J Neurosci 33(17):7463-7474.
%       [2] Raab DH (1962) Statistical facilitation of simple reaction
%           times. Trans NY Acad Sci 24(5):574-590.
%       [3] Miller J (1982) Divided attention: Evidence for coactivation
%           with redundant signals. Cogn Psychol 14(2):247-279.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 2-Apr-2019

% Decode input variable arguments
[q,per,lim,dep,test,area] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(4,2);
lims(1,:) = prctile(x,per);
lims(2,:) = prctile(y,per);
lims(3,:) = prctile(z,per);
lims(4,:) = prctile(xyz,per);

% Limit RTs to specified range
x = x(x>lims(1,1) & x<lims(1,2));
y = y(y>lims(2,1) & y<lims(2,2));
z = z(z>lims(3,1) & z<lims(3,2));
xyz = xyz(xyz>lims(4,1) & xyz<lims(4,2));

% Get min and max RT limits
if isempty(lim)
    lim = [min(lims(:)),max(lims(:))];
end

% Compute CDFs
fx = rt2cdf(x,q,lim);
fy = rt2cdf(y,q,lim);
fz = rt2cdf(z,q,lim);
fxyz = rt2cdf(xyz,q,lim);

% Compute Grice's bound
fmax = max([fx,fy,fz],[],2);

% Compute difference
if strcmpi(test,'ver')
    femp = fxyz-fmax;
elseif strcmpi(test,'hor')
    femp = fmax-fxyz;
end

% Compute empirical benefit
bemp = getauc(q,femp,area);

if nargout > 1
    
    % Compute race model
    if dep == 0 % Raab's Model
        fxy = fx+fy-fx.*fy;
        frace = fxy+fz-fxy.*fz;
    elseif dep == -1 % Miller's Bound
        frace = fx+fy+fz;
        frace(frace>1) = 1;
    end
    
    % Compute difference
    if strcmpi(test,'ver')
        fpred = frace-fmax;
    elseif strcmpi(test,'hor')
        fpred = fmax-frace;
    end
    
    % Compute predicted benefit
    bpred = getauc(q,fpred,area);
    
end

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
if any(strcmpi(varargin,'dep')) && ~isempty(varargin{find(strcmpi(varargin,'dep'))+1})
    dep = varargin{find(strcmpi(varargin,'dep'))+1};
    if dep~=0 && dep~=-1 && dep~=1
        error('DEP must be a scalar with a value of -1, 0 or 1.')
    end
else
    dep = 0; % default: assume statistical independence (Raab's Model)
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