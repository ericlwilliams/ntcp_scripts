function [sortQ,indxQ,indxorg]=EqualIntervals(numQ,numIntvl)
% numQ - the quantity vector which will be partitioned into intervals where each interval will have approximately same number of quantities.
% numIntvl - the number of intervals.
    
    [sortQ, indxorg] = sort(numQ(:));
    uniQ = unique(sortQ);
    h = hist(sortQ,uniQ);
    hcum = cumsum(h);
    intvsize = hcum(end)/numIntvl; % number of quantities in each interval
    
    indxQ = zeros(length(numQ(:)),1); % index for each element in numQ
    n=0; % the number of elements allocated
    for k=1:numIntvl % put the numQ into the intervals
        f = find(hcum>n+intvsize);
        if isempty(f)
            indxQ(n+1:end) = k;
            break;
        else
            if f(1)==1 % many elements have the same smallest quantity
                indxQ(n+1:n+hcum(f(1)))=k;
                n=n+hcum(f(1));
            elseif abs(hcum(f(1))-n - intvsize) > abs(hcum(f(1)-1)-n - intvsize)
                if hcum(f(1)-1) > n
                    indxQ(n+1:hcum(f(1)-1)) = k;
                    n = hcum(f(1)-1);
                else
                    indxQ(n+1:hcum(f(1))) = k;
                    n = hcum(f(1));
                end
            else
                indxQ(n+1:hcum(f(1))) = k;
                n = hcum(f(1));
            end
        end
        if n == hcum(end)
            break;
        else
            intvsize = (hcum(end)-n)/(numIntvl-k);
        end
    end
%     indxQ = indxQ(indx);