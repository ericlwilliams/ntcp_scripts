function GenerateDvhAtlas(analy_str,structure_str,tox_str,data_loc,atlas_loc)
%% Generates DVH atlases from institutional (or combined) meta data
%% control from runGenerateDvhAtlas: with GenerateDvhAtlas(analy_str,structure_str,tox_string);
% Examples
% analy_str = 'nfz','pld'
% structure_str = 'eso','ilung','heart','nfz','pbt','nfz','lungs'
% tox_str = 'rp','pultrox','esotox'
% data_loc = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta_comb.mat'
% atlas_loc = 'Z:/elw/MATLAB/meta_analy/atlases'

tic; %close all;
disp(['$$$$']);
disp(['Running GenerateDvhAtlas(',analy_str,',',structure_str,',',tox_str,',...',10,...
        data_loc,',...',10,...
        atlas_loc]);
disp(['$$$$']);

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%% Globals
dose_bin_step = 1; %Gy
vol_bin_step = 0.05; % steps of log10(v)

% TMP, run with no-arguments does pld analysis, data/arguments loaded below

if nargin==0,
    do_analysis='set analysis (or run with runGenerateDvhAtlas.m)'; %pld
else
    do_analysis=analy_str;
end
    
%% Analysis specific file loading
% Choose meta data to work with, one insitute (or combined) at a time

if isequal(do_analysis,'pld'),
%% Pooled analysis
atlas_loc = 'Z:/elw/MATLAB/meta_analy/atlases/';

% Combined
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta_comb.mat';
load(fn,'CGcomb'); CGdata = CGcomb; clear CGcomb; analy_str='rp_comb';

atlas_file_name = [atlas_loc,analy_str,'_aoc.xls'];

% MSKCC
%fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_fine_EUD_meta.mat';
%analy_str='rp_msk';

% NKI
%fn='Z:/elw/MATLAB/meta_analy/meta_data/NKI_fine_EUD_meta.mat';
%analy_str='rp_nki';

% RTOG
%fn='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_fine_EUD_meta.mat';
%analy_str='rp_rtog';

% UMich
%fn='Z:/elw/MATLAB/meta_analy/meta_data/UMich_fine_EUD_meta.mat';
%analy_str='rp_umich';

%load(fn,'CGobjs'); CGdata = CGobjs; clear CGobjs;

elseif isequal(do_analysis,'nfz') % specifics loaded from arguments
    nfzKeySet = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
    nfzValueSet = {'ilung','eso','heart','nfz','pbt','lungs'};
    nfzMapObj = containers.Map(nfzKeySet,nfzValueSet);
    a2b='Inf';

    % nfz meta data file to load
    file_name = ['NFZ_',structure_str,'_',tox_str,'_a2b',a2b,'_data.mat'];
    
    % set altas location/name
    atlas_file_name = [atlas_loc,analy_str,'_',...
                       nfzMapObj(structure_str),'_',...
                       tox_str,'_aoc.xls']; %atlas location and file name
    
    load(strcat(data_loc,file_name),'CGobj_org');
    CGdata = CGobj_org;
    clear CGobj_org;
    
elseif isequal(do_analysis,'eso')
    
     nfzKeySet = {'ESOPHAGUS'};
    nfzValueSet = {'eso'};
    nfzMapObj = containers.Map(nfzKeySet,nfzValueSet);
    a2b='Inf';

    % nfz meta data file to load
    file_name = ['NFZ_',structure_str,'_',tox_str,'_a2b',a2b,'_acute_data.mat'];
    
    % set altas location/name
    atlas_file_name = [atlas_loc,analy_str,'_',...
                       nfzMapObj(structure_str),'_',...
                       tox_str,'_aoc.xls']; %atlas location and file name
    
    load(strcat(data_loc,file_name),'CGobj_org');
    CGdata = CGobj_org;
    clear CGobj_org;
    
    
    
    
    
    
elseif isequal(do_analysis,'bpx')
    
    % set altas location/name
    atlas_file_name = [atlas_loc,analy_str,'_aoc.xls']; %atlas location and file name
    
    
    load(data_loc,'CGobj_current');
    CGdata = CGobj_current;
    clear CGobj_current;
    
else
      error(['!!! Analysis type: ',do_analysis,' unknown']); 
end

%% NFZ
% todo


%% Calculate generalized DVH dose bins
% TODO - check if physical (depends on inst.)

CGgrp = [CGdata.mGrp];
clear CGdata;

flgcensor = [CGgrp.mFlgCensor]; flgcomp = ~flgcensor;

%max dose
doses = cellfun(@(x) x(end), {CGgrp.mDoseBins_org});
dmax = max(doses);
dose_bins = ( 0 : dose_bin_step : dmax+dose_bin_step )';

% volume bins in equal steps of log10(v)

%if isempty([CGgrp.mVolCum]),
    for t=1:length(CGgrp)
        CGgrp(t) = CGgrp(t).fDiff2Cum();
    end
%end

mx_vols = cellfun(@(x) max(x), {CGgrp.mVolCum});
vmax = max(mx_vols);

%min non-zero
%min_vols = [CGgrp.mVolCum];
%min_vols = cellfun(@(x) min(x), {CGgrp.mVolCum});

% protect for unrealistically small min vols
min_vols = cellfun(@(x) min(x(x>0.1)), {CGgrp.mVolCum});
vmin = min(min_vols);
vmin = round(vmin/vol_bin_step)*vol_bin_step;%round to nearest bin step
vol_bins = [log10(vmin) : vol_bin_step : log10(vmax)+vol_bin_step];
vol_bins = [0, 10.^vol_bins];



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
                
                
                cur_vols(jj) = interp1( [FxGrp(jj).mDoseBins_org(f); FxGrp(jj).mDoseBins_org(f+1)], [FxGrp(jj).mVolCum(f); FxGrp(jj).mVolCum(f+1)], cur_dose );
                % normalize to percent maximum
                %cur_vols(jj) = 100*(cur_vols(jj)/max(FxGrp(jj).mVolCum));
            end
        end
        cur_vols(cur_vols==0) = -1; % exclude zero volume patients
        
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
    
    FxAtlas = cell(vol_dim+2,2*dose_dim+2);
    FxAtlas{1,1}='Vol (cc)';
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
    
    atlas_sheet_name = 'nfx = ';
    
    if fxs(i)<0, 
        atlas_sheet_name = [atlas_sheet_name,'all'];
    else
        atlas_sheet_name = [atlas_sheet_name,num2str(fxs(i))];
    end
    
warning('off','MATLAB:xlswrite:AddSheet');
%     testFile = xlsread(atlas_file_name,atlas_sheet_name);
%     if ~isempty(testFile)
%         xlswrite(atlas_file_name,zeros(size(testFile))*nan,atlas_sheet_name);
%     end
    disp(['Writing ',atlas_file_name,', sheet: ',atlas_sheet_name,'...']);
    xlswrite(atlas_file_name,FxAtlas,atlas_sheet_name,'A1');
    warning('on','MATLAB:xlswrite:AddSheet');
    

end

toc;
end