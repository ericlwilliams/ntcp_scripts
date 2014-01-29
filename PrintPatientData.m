%% INCOMPLETE
function PrintPatientData
tic; %close all;

%% load MSK/NKI data (outputs from Pt_MSK/NKI,

fn = 'Z:/elw/MATLAB/regions/data/EUD_regional_MSK_NKI.mat';
load(fn,'CGmsk','CGnki','CGcomb');

numreg = length(CGcomb)-2; %ignore left/right
n_msk_pts = zeros(1,numreg);
n_nki_pts = zeros(1,numreg);

for j=1:numreg,
    n_msk_pts(j) = CGmsk(j).mNumInGrp;
    n_nki_pts(j) = CGnki(j).numInGrp;
end

%% Print flags

printTx = true; % print dose delivered (treated)

if printTx,
    disp(['Dose delivered']);
    
    for i=1:numreg,
       mskGrp = CGmsk(1).mGrp;
       for j=1:n_msk_pts(i);
           disp([mskGrp(j).mID, ' -> ', mskGrp(j).mDoseTx]);
       end
    end
end



end
