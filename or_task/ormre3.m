function [mre,Fmre] = ormre3(x,y,z,xyz,p,varargin)
%ormre3 Multisensory response enhancement for a trisensory OR task.
%   MRE = ORMRE3(X,Y,Z,XYZ) returns the multisensory response enhancement
%   (MRE) for a trisensory OR task, quantified by the area between the CDFs
%   of the faster of the unisensory RT distributions X, Y and Z, and the
%   trisensory RT distribution XYZ, normalized by the former (Diederich &
%   Colonius, 2004). X, Y, Z and XYZ can have different lengths. This
%   function treats NaNs as missing values, and ignores them.
%
%   [...] = ORMRE3(...,P) uses the probabilities P to compute the CDF
%   quantiles. P is a vector of decimal values between 0 and 1 inclusive
%   (default=0.05:0.1:0.95).
%
%   [...,FMRE] = ORMRE3(...) returns the MRE at every CDF quantile. If the
%   mean or median method is chosen, then FMRE is returned empty.
%
%   [...] = ORMRE3(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([X,Y,Z,XYZ]),max([X,Y,Z,XYZ])])
%   'method'    a string specifying the method used to compute the MRE
%                   'cdf'       area between the CDFs (default)
%                   'mean'      mean RT values
%                   'median'    median RT values
%   'area'      a string specifying how to compute the area under the curve
%                   'all'       use all values (default)
%                   'pos'       use only positive values
%                   'neg'       use only negative values
%
%   See also ORMRE, ORBENEFIT3, ORGAIN3, ORMODEL3, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB
%           Package for Stochastic Modelling of Multisensory Reaction
%           Times (In prep).
%       [2] Diederich A, Colonius H (2004) Bimodal and trimodal
%           multisensory enhancement: Effects of stimulus onset and
%           intensity on reaction time. Percept Psychophys 66(8):1388-1404.
%       [3] Ulrich R, Miller J, Schroter H (2007) Testing the race model
%           inequality: An algorithm and computer programs. Behav Res
%           Methods 39(2):291-302.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 15-May-2019

% Decode input variable arguments
[lim,method,area] = decode_varargin(varargin);

% Set default values
if nargin < 4 || isempty(p)
    p = 0.05:0.1:0.95;
elseif ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
end

% Transpose row vectors
if isrow(x), x = x'; end
if isrow(y), y = y'; end
if isrow(z), z = z'; end
if isrow(xyz), xyz = xyz'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([x;y;z;xyz]),max([x;y;z;xyz])];
end

% Compute CDFs/CTs
if strcmpi(method,'cdf')
    Fx = rt2cfp(x,lim(2));
    Fy = rt2cfp(y,lim(2));
    Fz = rt2cfp(z,lim(2));
    Fxyz = rt2cfp(xyz,lim(2));
elseif strcmpi(method,'mean')
    Fx = nanmean(x);
    Fy = nanmean(y);
    Fz = nanmean(z);
    Fxyz = nanmean(xyz);
elseif strcmpi(method,'median')
    Fx = nanmedian(x);
    Fy = nanmedian(y);
    Fz = nanmedian(z);
    Fxyz = nanmedian(xyz);
end

% Compute CDF quantiles
if strcmpi(method,'cdf')
    Fx = cfp2q(Fx,p);
    Fy = cfp2q(Fy,p);
    Fz = cfp2q(Fz,p);
    Fxyz = cfp2q(Fxyz,p);
end

% Compute Grice's bound
Fbound = min([Fx,Fy,Fz],[],2);

% Compute multisensory response enhancement
Fmre = (Fbound-Fxyz)./Fbound*100;

% Compute MRE area
if strcmpi(method,'cdf')
    mre = getauc(p,Fmre,area);
else
    mre = Fmre;
    Fmre = [];
end

function [lim,method,area] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'lim')) && ~isempty(varargin{find(strcmpi(varargin,'lim'))+1})
    lim = varargin{find(strcmpi(varargin,'lim'))+1};
    if ~isnumeric(lim) || isscalar(lim) || any(isnan(lim)) || any(isinf(lim)) || any(lim<0) || lim(1)>=lim(2)
        error('LIM must be a 2-element vector of positive values.')
    end
else
    lim = []; % default: unspecified
end
if any(strcmpi(varargin,'method')) && ~isempty(varargin{find(strcmpi(varargin,'method'))+1})
    method = varargin{find(strcmpi(varargin,'method'))+1};
    if ~any(strcmpi(method,{'cdf','mean','median'}))
        error('Invalid value for argument METHOD. Valid values are: ''cdf'', ''mean'', ''median''.')
    end
else
    method = 'cdf'; % default: area between CDFs
end
if any(strcmpi(varargin,'area')) && ~isempty(varargin{find(strcmpi(varargin,'area'))+1})
    area = varargin{find(strcmpi(varargin,'area'))+1};
    if ~any(strcmpi(area,{'all','pos','neg'}))
        error('Invalid value for argument AREA. Valid values are: ''all'', ''pos'', ''neg''.')
    end
else
    area = 'all'; % default: use all values
end