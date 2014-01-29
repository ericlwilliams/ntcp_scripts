function pr=FisherExactTest(table2x2)
    if 1
        n1=sum(table2x2(1,:)); n2=sum(table2x2(2,:)); m1=sum(table2x2(:,1)); m2=sum(table2x2(:,2));
        if n1>n2 % exchange it so n1<n2
            table2x2=flipud(table2x2);
        end
        if m1>m2 % exchange so m1<m2
            table2x2=fliplr(table2x2);
        end
        n1=sum(table2x2(1,:)); n2=sum(table2x2(2,:)); m1=sum(table2x2(:,1)); m2=sum(table2x2(:,2)); n=m1+m2;
    else
        table2x2=sortrows(table2x2); table2x2= sortrows(table2x2'); table2x2=table2x2'; % put the smallest value in the first row and first column
        table2x2(1,2)=n1; table2x2(2,1)=m1; table2x2(2,2)=n2-m1;
    end
    pr=zeros(min(n1,m1)+1,1);
    % the first probability
    d1 = n2-m1+1:n2; d2 = m2+1:n;
    d1 = d1./d2; pr(1)=prod(d1);
    % the rest probability
    for k=1:min(n1,m1)
        pr(k+1) = pr(k) * ( (table2x2(1,2)-k+1) * (table2x2(2,1)-k+1) / ( k * (table2x2(2,2)+k)));
    end
