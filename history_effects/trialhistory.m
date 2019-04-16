function [x,y,xy,labels] = trialhistory(data,cond,varargin)
%trialhistory Separate trials based on their trial history.
%   [X,Y,XY] = TRIALHISTORY(DATA,COND) returns the trials from DATA that
%   were preceded by a different condition (i.e., switch trials). DATA
%   should be a vector of mixed unisensory and multisensory reaction times,
%   arranged in the order in which they were collected during testing. COND
%   should be a vector of integer values indicating the condition of each
%   trial, where 1 and 2 indicate unisensory trials, and 3 indicates
%   multisensory trials.
%
%   [...,LABELS] = TRIALHISTORY(...) returns a vector of labels indicating
%   the condition of both the previous and current trial.
%
%   Trials are labelled as follows:
%
%       Condition   Label   Type
%       X->X        1       repeat
%       Y->X        2       switch
%       XY->X       3       repeat*
%       X->Y        4       switch
%       Y->Y        5       repeat
%       XY->Y       6       repeat*
%       X->XY       7       switch**
%       Y->XY       8       switch**
%       XY->XY      9       repeat
%
%   *not a pure repeat
%   **not a pure switch
%
%   [...] = TRIALHISTORY(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'type'      a string specifying the type of trial history
%                   'switch'    switch trials (default)
%                   'repeat'    repeat trials
%                   'x'         trials preceded by condition X
%                   'y'         trials preceded by condition Y
%                   'xy'        trials preceded by condition XY
%   'uni'       a scalar specifying whether to include unisensory trials
%               preceded by multisensory trials (*): pass in 0 to exclude
%               (default) or 1 to include
%   'multi'     a scalar specifying whether to include multisensory trials
%               preceded by both unisensory conditions (**): pass in 0 to
%               include both (default), 1 to include only the first
%               unisensory condition, or 2 to include only the second
%               unisensory condition
%
%   See also SWITCHCOST, BIASMODEL, BIASGAIN, BIASBENEFIT, F1SCORE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 16-Apr-2019

% Decode input variable arguments
[type,multi,uni] = decode_varargin(varargin);

% Get labels for each trial
labels = zeros(length(data),1);
for i = 1:length(data)
    if i == 1
        labels(i) = NaN;
    elseif cond(i-1) == 1 && cond(i) == 1 % 1: X->X
        labels(i) = 1;
    elseif cond(i-1) == 2 && cond(i) == 1 % 2: Y->X
        labels(i) = 2;
    elseif cond(i-1) == 3 && cond(i) == 1 % 3: XY->X
        labels(i) = 3;
    elseif cond(i-1) == 1 && cond(i) == 2 % 4: X->Y
        labels(i) = 4;
    elseif cond(i-1) == 2 && cond(i) == 2 % 5: Y->Y
        labels(i) = 5;
    elseif cond(i-1) == 3 && cond(i) == 2 % 6: XY->Y
        labels(i) = 6;
    elseif cond(i-1) == 1 && cond(i) == 3 % 7: X->XY
        labels(i) = 7;
    elseif cond(i-1) == 2 && cond(i) == 3 % 8: Y->XY
        labels(i) = 8;
    elseif cond(i-1) == 3 && cond(i) == 3 % 9: XY->XY
        labels(i) = 9;
    end
end

% Separate trials by previous condition
if strcmpi(type,'repeat')
    x = data(labels==5);
    y = data(labels==9);
    xy = data(labels==1);
elseif strcmpi(type,'switch')
    if uni == 0 % Exclude **conditions
        x = data(labels==6);
        y = data(labels==8);
    elseif uni == 1 % Include **conditions
        x = data(labels==4 | labels==6);
        y = data(labels==7 | labels==8);
    end
    if multi == 0 % Include both *conditions
        xy = data(labels==2 | labels==3);
    elseif multi == 1 % Include first *condition
        xy = data(labels==2);
    elseif multi == 2 % Include second *condition
        xy = data(labels==3);
    end
elseif strcmpi(type,'x')
    x = data(labels==1);
    y = data(labels==4);
    xy = data(labels==7);
elseif strcmpi(type,'y')
    x = data(labels==2);
    y = data(labels==5);
    xy = data(labels==8);
elseif strcmpi(type,'xy')
    x = data(labels==3);
    y = data(labels==6);
    xy = data(labels==9);
end

function [type,multi,uni] = decode_varargin(varargin)
%decode_varargin Decode input variable arguments.
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'type')) && ~isempty(varargin{find(strcmpi(varargin,'type'))+1})
    type = varargin{find(strcmpi(varargin,'type'))+1};
    if ~any(strcmpi(type,{'ver','hor'}))
        error('Invalid value for argument TYPE. Valid values are: ''switch'', ''repeat'', ''x'', ''y'', ''xy''.')
    end
else
    type = 'switch'; % default: switch trials
end
if any(strcmpi(varargin,'uni')) && ~isempty(varargin{find(strcmpi(varargin,'uni'))+1})
    uni = varargin{find(strcmpi(varargin,'uni'))+1};
    if uni~=0 && uni~=1
        error('UNI must be a scalar with a value of 0 or 1.')
    end
else
    uni = 0; % default: exclude *conditions
end
if any(strcmpi(varargin,'multi')) && ~isempty(varargin{find(strcmpi(varargin,'multi'))+1})
    multi = varargin{find(strcmpi(varargin,'multi'))+1};
    if multi~=0 && multi~=1 && multi~=2
        error('MULTI must be a scalar with a value of 0, 1 or 2.')
    end
else
    multi = 0; % default: include both **conditions
end