function [fx,fy,fxy,frace,fdiff] = racemodel(x,y,xy,varargin)
%racemodel Unisensory, multisensory and race model reaction time CDFs.
%   [FX,FY,FXY] = RACEMODEL(X,Y,XY) returns the cumulative distribution
%   functions (CDFs) for the unisensory RT distributions X and Y, and the 
%   multisensory RT distribution XY at 20 linearly-spaced quantiles between 
%   0.05 and 1.
% 
%   [...,FRACE] = RACEMODEL(...) returns the CDF of the race model based on 
%   the unisensory RT distributions X and Y. The race model is computed 
%   using probability summation (Raab, 1962), which assumes statistical 
%   independence between X and Y. For valid estimates of FRACE, the stimuli 
%   used to generate X, Y and XY should be randomly interleaved in order 
%   for the assumption of context invariance to hold (Luce, 1986).
% 
%   [...,FDIFF] = RACEMODEL(...) returns the difference between FXY and 
%   FRACE to test whether XY exceeded statistical facilitation and thus 
%   "violated" the race model (Miller, 1982).
% 
%   [...] = RACEMODEL(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'per'       a vector specifying the range of quantiles used to compute
%               the CDFs (default=[0.05:0.05:1])
%   'range'     a 2-element vector specifying the lower and upper
%               percentiles of RTs to consider (default=[0,100])
%   'dep'       a scalar specifying whether statistical dependence between
%               X and Y is assumed: pass in 0 to assume independence (Raab, 
%               1962; default) or 1 to assume dependence (Miller, 1982)
%   'test'      a string specifying how to test the race model
%                   'ver'       vertical test (default)
%                   'hor'       horizontal test 
% 
%   See also RSEGAIN, RSEBENEFIT, PERMTEST_T, EFFECTSIZE_D.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Raab DH (1962) Statistical facilitation of simple reaction 
%           times. Trans NY Acad Sci, 24(5):574-590.
%       [2] Luce RD (1986) Response times: Their role in inferring mental
%           organization. New York, NY: Oxford University Press.
%       [3] Miller J (1982) Divided attention: Evidence for coactivation 
%           with redundant signals. Cogn Psychol, 14(2):247-279.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 5-Feb-2019

% Decode input variable arguments
[per,range,dep,test] = decode_varargin(varargin);

% Get RT range for each condition
lims = zeros(3,2);
lims(1,:) = prctile(x,range);
lims(2,:) = prctile(y,range);
lims(3,:) = prctile(xy,range);

% Limit RTs to specified range 
x = x(x>lims(1,1) & x<lims(1,2));
y = y(y>lims(2,1) & y<lims(2,2));
xy = xy(xy>lims(3,1) & xy<lims(3,2));

% Compute linearly-spaced quantiles between min/max RTs
qntls = linspace(min(lims(:)),max(lims(:)),length(per));

% Compute cumulative distribution functions
fx = zeros(length(per),1); nx = length(x);
fy = zeros(length(per),1); ny = length(y);
fxy = zeros(length(per),1); nxy = length(xy);
for i = 1:length(per)
    fx(i) = sum(x<=qntls(i))/nx;
    fy(i) = sum(y<=qntls(i))/ny;
    fxy(i) = sum(xy<=qntls(i))/nxy;
end

% Compute race model
if nargout > 3
    if dep == 0
        frace = fx+fy-fx.*fy;
    elseif dep == 1
        frace = fx+fy; 
        frace(frace>1) = 1;
    end
end

% Compute race model violation
if nargout > 4
    fdiff = fxy-frace;
end

function [per,range,dep,test] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'per')) && ~isempty(varargin{find(strcmpi(varargin,'per'))+1})
    per = varargin{find(strcmpi(varargin,'per'))+1};
    if ~isnumeric(per) || isscalar(per) || any(isnan(per)) || any(isinf(per)) || any(per<0) || any(per>1) || per(1)>=per(2)
        error('PER must be a vector with values between 0 and 1.')
    end
else
    per = 0.05:0.05:1; % default: 0.05 to 1 in 0.05 increments
end
if any(strcmpi(varargin,'range')) && ~isempty(varargin{find(strcmpi(varargin,'range'))+1})
    range = varargin{find(strcmpi(varargin,'range'))+1};
    if ~isnumeric(range) || isscalar(range) || any(isnan(range)) || any(isinf(range)) || any(range<0) || any(range>100) || range(1)>=range(2)
        error('RANGE must be a 2-element vector with values between 0 and 100.')
    end
else
    range = [0,100]; % default: all RTs
end
if any(strcmpi(varargin,'dep')) && ~isempty(varargin{find(strcmpi(varargin,'dep'))+1})
    dep = varargin{find(strcmpi(varargin,'dep'))+1};
    if any(dep~=0) && any(dep~=1)
        error('DEP must be a scalar with a value of 0 (independence) or 1 (dependence).')
    end
else
    dep = 0; % default: assume statistical independence
end
if any(strcmpi(varargin,'test')) && ~isempty(varargin{find(strcmpi(varargin,'test'))+1})
    test = varargin{find(strcmpi(varargin,'test'))+1};
    if ~any(strcmpi(test,{'ver','hor'}))
        error('Invalid value for argument TEST. Valid values are: ''ver'', ''hor''.')
    end
else
    test = 'ver'; % default: vertical test
end