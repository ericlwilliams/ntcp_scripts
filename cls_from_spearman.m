function [cl68_high,cl68_low,cl95_high,cl95_low]=cls_from_spearman(spear_r,ndf)

    z=0.5*log((1+spear_r)/(1-spear_r));

    z68_low=(spear_r-(1/sqrt(ndf-3)));
    z68_high=(spear_r+(1/sqrt(ndf-3)));

    z95_low=(spear_r-(1.96/sqrt(ndf-3)));
    z95_high=(spear_r+(1.96/sqrt(ndf-3)));

    cl68_low = (exp(2*z68_low)-1)/(exp(2*z68_low)+1);
    cl68_high = (exp(2*z68_high)-1)/(exp(2*z68_high)+1);

    cl95_low = (exp(2*z95_low)-1)/(exp(2*z95_low)+1);
    cl95_high = (exp(2*z95_high)-1)/(exp(2*z95_high)+1);

end
