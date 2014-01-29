function chi2=Chi2From2by2Table(table2x2)
    n1=sum(table2x2(1,:)); n2=sum(table2x2(2,:)); m1=sum(table2x2(:,1)); m2=sum(table2x2(:,2)); n=m1+m2;
    chi2 = (n-1) * (table2x2(1,1)*table2x2(2,2) - table2x2(2,1)*table2x2(1,2))^2;
    chi2 = chi2 / (n1*n2*m1*m2);
