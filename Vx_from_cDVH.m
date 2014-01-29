function vx=Vx_from_cDVH(cDVH,x_dose)
% Calculate geud from given absolute dDVH
% converte absolute DVH to relative DVH

tot_vol = cDVH(1,2);
relDVH = cDVH;
relDVH(:,2) = relDVH(:,2)./tot_vol;% abs -> rel
relDVH(:,1) = relDVH(:,1)./100.;% cGy->Gy
vx=0;
for i=2:2:length(cDVH(:,1)),
    if x_dose<=relDVH(i,1),
     vx=relDVH(i,2);
     break;
    end
end
vx=vx*100;


end