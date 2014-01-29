function geud=gEUD_from_absDVH(absDVH,a)
% Calculate geud from given absolute dDVH
% converte absolute DVH to relative DVH

tot_vol = sum(absDVH(2:2:end,2));
%message = sprintf('tot_vol = %g cc',tot_vol);
%disp(message);disp(' ');

relDVH = absDVH;
relDVH(:,2) = relDVH(:,2)./tot_vol;
relDVH(:,1) = relDVH(:,1)./100.;
%assignin('base','relDVH',relDVH);
geud = (sum((relDVH(2:2:end,1).^a).*relDVH(2:2:end,2)))^(1/a);


end