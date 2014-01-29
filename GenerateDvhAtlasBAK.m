function GenerateDvhAtlas
%% Generates DVH atlases from institutional (or combined) meta data
tic; %close all;

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%% Globals
dose_bin_step = 1; %Gy
vol_bin_step = 5; % percentage

%% Analysis specific file loading
% Choose meta data to work with, one insitute (or combined) at a time

% Pooled analysis
%fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_me
%ta.mat';

% Combined
%load(fn,'CGcomb'); CGdata = CGcomb; clear CGcomb;
% MSKCC
%load(fn,'CGmsk'); CGdata = CGmsk; clear CGmsk;
% NKI
%load(fn,'CGnki'); CGdata = CGnki; clear CGnki;
% RTOG
%load(fn,'CGrtog'); CGdata = CGrtog; clear CGrtog;
% UMich
%load(fn,'CGumich'); CGdata = CGumich; clear CGumich;

%% Testing ONLY, REMOVE
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_fine_EUD_meta.mat'; %% Testing only!
load(fn,'CGobjs'); CGdata = CGobjs; clear CGobjs;
atlas_loc = 'Z:/elw/MATLAB/meta_analy/atlases/';
analy_str ='rp_mskcc_aoc';

%% NFZ
% todo


%% Calculate generalized DVH dose bins
% TODO - check if physical (depends on inst.)

CGgrp = [CGdata.mGrp];
clear CGdata;

flgcensor = [CGgrp.mFlgCensor]; flgcomp = ~flgcensor;

% only interested in highest dose bin from complication group
doses = cellfun(@(x) x(end), {CGgrp(flgcomp).mDoseBins_org});
dmax = max(doses);
dose_bins = ( 0 : dose_bin_step : dmax+dose_bin_step )';

vol_bins = [0:vol_bin_step:100];


% will write atlas for each fraction bin and for combined
grp_flgs = [CGgrp.mFxNum];

%unique fractionations
fxs = unique(grp_flgs);
% get number of entries for each unique fx number
fx_counts = histc(grp_flgs,fxs);

% add combined fractionation schemes
fxs = [-1 fxs];
fx_counts = [sum(fx_counts) fx_counts];

%% loop over each fractionation scheme (last entry is combined schemes)
for i=1:length(fxs),
    cur_fx = fxs(i);
    cur_fx_pts = fx_counts(i);   
    fx_flg = true(length(CGgrp),1)';
    if ~(cur_fx<0), % cur_fx = -1; combines all fxs
        fx_flg = logical([CGgrp.mFxNum]==cur_fx);
    end
    FxGrp = CGgrp(fx_flg);
    FxComp = flgcomp(fx_flg);
    
    %% atlas generation
    
    AtlasTotal = zeros( length(dose_bins), length(vol_bins));
    AtlasComp = zeros( length(dose_bins), length(vol_bins));
    
    %% Rebin patient DVHs given dose_bins
    
    
    for ii=1:length(dose_bins),
        
        cur_dose = dose_bins(ii);
        cur_vols = zeros(cur_fx_pts,1); % patient volumes for current dose bin
        
        for jj=1:cur_fx_pts,
            
           f = find( FxGrp(jj).mDoseBins_org <= cur_dose ); f=f(end); % f won't be empty by definition
            
            if FxGrp(jj).mDoseBins_org(f) < FxGrp(jj).mDoseBins_org(end) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero
                
                if isempty(FxGrp(jj).mVolCum),
                    FxGrp(jj) = FxGrp(jj).fDiff2Cum();
                end
                
                cur_vols(jj) = interp1( [FxGrp(jj).mDoseBins_org(f); FxGrp(jj).mDoseBins_org(f+1)], [FxGrp(jj).mVolCum(f); FxGrp(jj).mVolCum(f+1)], cur_dose );
                % normalize to percent maximum
                cur_vols(jj) = 100*(cur_vols(jj)/max(FxGrp(jj).mVolCum));
            end
        end
        %cur_vols(cur_vols==0) = -1; % exclude zero volume patients
        
        for kk=1:length(vol_bins)
            vol_f = find( cur_vols >= vol_bins(kk)); % patient at the grid point
            vol_g = find( FxComp(vol_f) );
            
            AtlasTotal(ii,kk) = length(vol_f);
            AtlasComp(ii,kk) = length(vol_g);
            
        end
    end
    %% Save current atlas to file
    % following convention from Fan in Dose-Volume Parameters Predit for
    % the Development of Chest Wall Pain After Stereotactic Body Radiation
    % for Lung Cancer, Mutter et al.
    
    % to print: flipud(AtlasTotal')
    % convert to cell array
    vol_dim = length(vol_bins);
    dose_dim = length(dose_bins);    
    
    % TMP - using rel vol, change to abs
    FxAtlas = cell(vol_dim+2,2*dose_dim+2);
    FxAtlas{1,1}='Vol (%)';
    FxAtlas{1,2}='log10(vol)\Dose (Gy)';
    FxAtlas(2:end-1,1)=mat2cell(vol_bins',ones(vol_dim,1),1);
    FxAtlas(2:end-1,2)=mat2cell(log10(vol_bins)',ones(vol_dim,1),1);
    FxAtlas{2,2}=[]; % remove if not using log10(v)
    FxAtlas(1,3:2:end) = mat2cell(dose_bins',1,ones(dose_dim,1));
    
    % add data
    FxAtlas(2:end-1,3:2:end) = mat2cell(AtlasComp',ones(vol_dim,1),ones(dose_dim,1));
    FxAtlas(2:end-1,4:2:end) = mat2cell(AtlasTotal',ones(vol_dim,1),ones(dose_dim,1));
    
    FxAtlas(end,3:2:end) = {'nc(D(i),V(j))'};
    FxAtlas(end,4:2:end) = {'np(D(i),V(j))'};
    
    atlas_file_name = [atlas_loc,analy_str,'.xls'];
    
    atlas_sheet_name = 'nfx = ';
    
    if fxs(i)<0,
        atlas_sheet_name = [atlas_sheet_name,'all'];
    else
        atlas_sheet_name = [atlas_sheet_name,num2str(fxs(i))];
    end
    
    %warning('off','MATLAB:xlswrite:AddSheet');
    xlswrite(atlas_file_name,FxAtlas,atlas_sheet_name,'A1');
    %warning('on','MATLAB:xlswrite:AddSheet');
    
    
end

%             pt = CGgrp;
%             num = length(pt); % number of patients

%             % Vx
%             vols = zeros( num, 1 );
%             for ii = 1:length(dose_bins)
%                 % volume of each patient at current dose (x)
%                 vols(:) = 0;
%                 for jj = 1:num
%                     cur_vol = pt(jj).fVolAtDose( dose_bins(ii) );
%                     vols(jj) = cur_vol/max(pt(jj).mVolCum);
%                     vols(jj) = vols(jj).*100;
%                 end
%                 vols(vols==0) = -1; % exclude zero volume patients
%
%                 % matrix at each (Di, Vj, Tk)
%                 for kk = 1:length(vol_bins) % for each volume point under dose x
%                     %f = find( vols >= CGdata.mBinsVol(jj) ); % patient at the grid point
%                 %f = find( vols >= (CGdata.mBinsVol(jj)/max(CGdata.mBinsVol)) ); % patient at the grid point
%
%                 f = find( vols >= vol_bins(kk)); % patient at the grid
%                 point
%                     g = find( flgcomp(f) );
%
%                     % for each time point
%                         % total patients
%                         CGdata.mAtlasTotal_DVH(ii,kk) = length(f);
%                         % patients with complications
%                         CGdata.mAtlasComp_DVH(ii,kk) = length(g);
%                 end
%             end
%


toc;
end