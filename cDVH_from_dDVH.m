function cDVH=cDVH_from_dDVH(dDVH)
% cDVH_from_dDVH returns culmaltive DVH as calculated from differential
clear cDVH;
for i=1:length(dDVH(:,1)),
    cDVH(i,1) = dDVH(i,1);
    cDVH(i,2) = sum(dDVH(i:2:end,2));
end
end