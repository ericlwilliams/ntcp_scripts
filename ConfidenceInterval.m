function [pLoc, CI] = ConfidenceInterval(xrange,loglikely,pct)
% find the boundary of confidence interval of the peak
% xrange - the parameter coordinates
% loglikely - the statistic log-likelihood curve of the xrange
% pct - the confidence level
% pLoc - the peak location
% CI - the confidence interval, interpolated linearly if not fall on
% exactly index.

% prepare
    CI = zeros(2,1);
    xrange = xrange(:);
    loglikely = loglikely(:);

% find the maxium point position in xrange
    [~,Loc] = max(loglikely);
    pLoc = xrange(Loc);

% the left end point
    f = find(loglikely(1:Loc)<=pct);
    if isempty(f)
        CI(1) = -inf; % left end is not defined
    else
        CI(1) = interp1(loglikely(f(end):f(end)+1), xrange(f(end):f(end)+1), pct);
    end
% the right end point
    f = find(loglikely(Loc:end)<=pct);
    if isempty(f)
        CI(2) = inf;
    else
        CI(2) = interp1(loglikely(Loc+f(1)-2:Loc+f(1)-1), xrange(Loc+f(1)-2:Loc+f(1)-1), pct);
    end
end