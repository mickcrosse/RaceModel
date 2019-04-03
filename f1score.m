function [f1,hitrate,errorrate] = f1score(x,rtmin,rtmax,nresp,ntarget)
%f1score F1 score of a test's detection accuracy.
%   F1 = F1SCORE(X,RTMIN,RTMAX,NRESP,NSTIM) returns the F1 score of a
%   test's detection accuracy using the RTs contained in X. RTs below RTMIN
%   are considered false alarms and RTs above RTMAX are considered misses.
%   The total number of responses NRESP should include all false alarms and
%   double-presses. Misses, whereby no response was registered, can be
%   included in X as missing values (NaNs). The F1 score is computed as the
%   harmonic mean of precision and recall (Van Rijsbergen, 1979).
%
%   [...,HITRATE] = F1SCORE(...) returns the hit rate based on the ratio of
%   hits to targets without consideration of false alarms.
%
%   [...,ERRORRATE] = F1SCORE(...) returns the error rate based on the
%   ratio of false alarms to targets.
%
%   [...] = F1SCORE(...,NTARGET) computes the F1 score by specifying the
%   exact number of targets NTARGET. This is useful if NTARGET does not
%   equal the number of elements in X because misses were not included in
%   X as NaNs or if RTs were removed due to outlier correction procedures.
%
%   See also RACEMODEL, RSEGAIN, RSEBENEFIT, TPERMTEST, EFFECTSIZE.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   References:
%       [1] Van Rijsbergen CJ (1979) Information Retrieval: Butterworth-
%           Heinemann.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 3-Apr-2019

% Set default values
if nargin < 5 || isempty(ntarget)
    ntarget = length(x);
end

% Find number of hits and misses
noRT = isnan(x) | x>rtmax;
misses = sum(noRT)+ntarget-length(x);
hits = ntarget-misses;

% Compute number of false alarms
falarms = sum(nresp)-hits+sum(x<rtmin);

% Compute precision and recall
precis = hits./(hits+falarms);
recall = hits./(hits+misses);

% Compute F1 score
f1 = 2*(precis.*recall)./(precis+recall);

% Compute hit rate
if nargout > 1
    hitrate = hits/ntarget;
end

% Compute false alarm rate
if nargout > 2
    errorrate = falarms/ntarget;
end