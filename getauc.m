function auc = getauc(x,y,p)
%getauc Get area under the curve.
%   AUC = GETAUC(X,Y,P) returns the area under the curve Y with respect to
%   X based on the portion of the curve P using trapezoidal numerical
%   integration. Valid values for argument P are 'all' (entire portion,
%   default), 'pos' (positive portion), and 'neg' (negative portion). X and
%   Y must be vectors of equal length.
%
%   See also RSEGAIN, RSEBENEFIT, SWITCHCOST, RT2CDF.
%
%   RaceModel https://github.com/mickcrosse/RaceModel

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Apr 2017; Last Revision: 3-Apr-2019

% Set default values
if nargin < 3 || isempty(p)
    p = 'all';
end
if nargin < 2 || isempty(y)
    y = x;
    x = 1/length(y):1/length(y):1;
end

% Convert row vectors to column vectors
if size(x,2)>1, x = x'; end
if size(y,2)>1, y = y'; end

% Compute x-intercepts
if ~strcmpi(p,'all')
    
    % Invert signal to get negative area
    if strcmpi(p,'neg')
        y = -y;
    end
    
    % Compute x-intercepts
    z = x;
    neg = false(length(y),1);
    for i = 1:length(y)-1
        if y(i)<0 && y(i+1)>0 && neg(i)==0
            z(i) = x(i)-y(i)/(y(i+1)-y(i))*diff(x(i:i+1));
            neg(i) = true;
        elseif y(i)>0 && y(i+1)<0
            z(i+1) = x(i)-y(i)/(y(i+1)-y(i))*diff(x(i:i+1));
            neg(i+1) = true;
        elseif y(i)<0 && y(i+1)>0 && neg(i)==1
            z = [z(1:i);x(i)-y(i)/(y(i+1)-y(i))*diff(x(i:i+1));z(i+1:end)];
            x = [x(1:i);x(i)-y(i)/(y(i+1)-y(i))*diff(x(i:i+1));x(i+1:end)];
            y = [y(1:i);0;y(i+1:end)];
            neg = [neg(1:i);true;neg(i+1:end)];
        end
    end
    x = z;
    
    % Replace negative values with zeros
    neg(y<=0) = true;
    y(neg) = 0;
    
    % Invert signal again for negative area
    if strcmpi(p,'neg')
        y = -y;
    end
    
end

% Compute AUC using trapezoidal method
auc = trapz(x,y);