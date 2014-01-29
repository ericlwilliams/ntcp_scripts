function dx=Dx_from_cDVH(cDVH,v_frac)
% Calculate geud from given absolute dDVH
% converte absolute DVH to relative DVH
v_frac=v_frac/100;
tot_vol = cDVH(1,2);
relDVH = cDVH;
relDVH(:,2) = relDVH(:,2)./tot_vol;
relDVH(:,1) = relDVH(:,1)./100.;
dx=0;
for i=2:2:length(cDVH(:,2)),
    if v_frac>=relDVH(i,2),
     dx=relDVH(i,1);
     break;
    end
end

end