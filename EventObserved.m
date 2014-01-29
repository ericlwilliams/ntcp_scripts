function [medianVal,meanVal,stdMeanVal,binLow,binHigh,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flgCensor,val,numIntervals)
    [sortQ,indxQ,indxorg] = EqualIntervals(val,numIntervals); % group patients according to their gEUDs
    
    medianVal = zeros(numIntervals,1);
    meanVal = zeros(numIntervals,1);
    stdMeanVal = zeros(numIntervals,1);
    binLow = zeros(numIntervals,1);
    binHigh = zeros(numIntervals,1);
    numComp = zeros(numIntervals,1);
    numTotal= zeros(numIntervals,1);
%             stdprob = zeros(numIntervals,1);
%             BetaInv84 = zeros(numIntervals,1);
%             BetaInv16 = zeros(numIntervals,1);
    flg = ~flgCensor(indxorg);
    for m = 1 : numIntervals
        binLow(m) = min(sortQ(indxQ==m));
        binHigh(m) = max(sortQ(indxQ==m));
        medianVal(m) = median(sortQ(indxQ==m));
        meanVal(m) = mean(sortQ(indxQ==m));
        stdMeanVal(m) = std(sortQ(indxQ==m));
        numComp(m) = sum(flg(indxQ==m));
        numTotal(m) = length(find(indxQ==m));
    end
%                 stdprob(m) = sqrt(numTotal*prob(m)*(1-prob(m)))/numTotal;
    BetaInv84 = betainv( .84, numComp+1, numTotal - numComp + 1 );
    BetaInv16 = betainv( .16, numComp+1, numTotal - numComp + 1 );
end
