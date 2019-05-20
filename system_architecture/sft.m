function [mic,sic,t] = sft(hh,ll,hl,lh,p,varargin)
%sft Systems factorial technology.
%   MIC = SFT(HH,LL,HL,LH) returns the mean interaction contrast (MIC)
%   based on the mean values of the bisensory conditions HH, LL, HL and LH, 
%   where H indicates high salience and L indicates low salience on the 
%   corresponding sensory channels. HH, LL, HL and LH can have different
%   lengths. This function treats NaNs as missing values, and ignores them.
%
%   [...,SIC] = SFT(...) returns the survivor interaction contrast (SIC) 
%   based on the cumulative distributions functions (CDFs) of the bisensory 
%   conditions HH, LL, HL and LH.
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
%   [...,T] = SFT(...) returns the time intervals used to compute the CDFs.
%
%   [...] = SFT(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values. Valid parameters are the following:
%
%   Parameter   Value
%   'lim'       a 2-element vector specifying the lower and upper RT limits
%               for computing CDFs: it is recommended to leave this
%               unspecified unless comparing directly with other conditions
%               (default=[min([HH,LL,HL,LH]),max([HH,LL,HL,LH])])
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
%   Apr 2017; Last Revision: 3-May-2019

% Decode input variable arguments
[lim,verbose] = decode_varargin(varargin);

% Set default values
if nargin < 5 || isempty(p)
    p = 0.05:0.1:0.95;
elseif ~isnumeric(p) || isscalar(p) || any(p<0|p>1)
    error('P must be a vector of values between 0 and 1.')
end

% Transpose row vectors
if isrow(hh), hh = hh'; end
if isrow(ll), ll = ll'; end
if isrow(hl), hl = hl'; end
if isrow(lh), lh = lh'; end

% Get min and max CDF limits
if isempty(lim)
    lim = [min([hh;ll;hl;lh]),max([hh;ll;hl;lh])];
end

% Compute CDFs
Fxhh = rt2cdf(hh,p,lim);
Fxhl = rt2cdf(hl,p,lim);
Fxlh = rt2cdf(lh,p,lim);
[Fxll,t] = rt2cdf(ll,p,lim);

% Compute mean interaction contrast
mic = nanmean(hh)+nanmean(ll)-nanmean(hl)-nanmean(lh);

% Compute survivor interaction contrast
sic = Fxhl+Fxlh-Fxhh-Fxll;

% Display reference table
if verbose
    Model = {'Serial, OR';'Serial, AND';'Parallel, OR';'Parallel, AND';'Coactive'};
    MIC = {'0';'0';'>0';'<0';'>0'};
    SIC = {'0';'-to+';'>0';'<0';'-to+'};
    reftable = table(Model,MIC,SIC);
    disp(reftable)
end

function [lim,verbose] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'verbose')) && ~isempty(varargin{find(strcmpi(varargin,'verbose'))+1})
    verbose = varargin{find(strcmpi(varargin,'verbose'))+1};
    if verbose~=0 && verbose~=1
        error('VERBOSE must be a scalar with a value of 0 or 1.')
    end
else
    verbose = 0; % default: display table
end