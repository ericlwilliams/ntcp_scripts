function GenerateDvhActuarialAtlas(analy_str,structure_str,tox_str,data_loc,atlas_loc)
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
atlas_loc = 'Z:/elw/MATLAB/meta_analy/atlases/latest/';

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
%     atlas_file_name = [atlas_loc,analy_str,'_',...
%                        nfzMapObj(structure_str),'_',...
%                        tox_str,'_aoc.xls']; %atlas location and file name
     atlas_file_name = [atlas_loc,'clung_',...
                        nfzMapObj(structure_str),'_',...   
                       tox_str,'_'];
                       
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


%%% 
f2 = ~cellfun('isempty',{CGgrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{CGgrp.mDateLastFollowup}); % patients with no last follow up date
compdate = inf(CGdata.mNumInGrp,1);
lastfollowup = inf(CGdata.mNumInGrp,1);
compdate(f2) = ([CGgrp(f2).mDateComp] - [CGgrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([CGgrp(f3).mDateLastFollowup] - [CGgrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );

clear CGdata;

%time_points = [-1 3 6 9 12 15 18 21 50];
time_points = [3 6 9 12 15 18 21 50 round(max(compdate+0.5))];



%loop over all fractions
% will write atlas for each fraction bin and for combined
grp_flgs = [CGgrp.mFxNum];

%unique fractionations
fxs = unique(grp_flgs);

% add combined fractionation schemes
fxs = [-1 fxs];

% loop over time points
for k=1:length(time_points)
    
    cur_time_point = time_points(k);
    if cur_time_point<0
       cur_atlas_sheet_name=['t=all']; 
       onsheet = true(length(compdate),1);
    else
        cur_atlas_sheet_name=['t=',num2str(cur_time_point),'m'];
        %atrisk = [compdate>=cur_time_point];
        
        % for a time t, want patients who either had a comp or were
        % censored BEFORE time t

        onsheet = [compdate<cur_time_point];
    end
    
    %cur_atlas_file_name = [atlas_file_name,'t',num2str(time_points(k)),'m'];
    
    

curCGgrp = CGgrp(onsheet);

flgcensor = [curCGgrp.mFlgCensor]; flgcomp = ~flgcensor;


%max dose
doses = cellfun(@(x) x(end), {curCGgrp.mDoseBins_org});
dmax = max(doses);
dose_bins = ( 0 : dose_bin_step : dmax+dose_bin_step )';

% volume bins in equal steps of log10(v)

%if isempty([curCGgrp.mVolCum]),
    for t=1:length(curCGgrp)
        curCGgrp(t) = curCGgrp(t).fDiff2Cum();
    end
%end

mx_vols = cellfun(@(x) max(x), {curCGgrp.mVolCum});
vmax = max(mx_vols);

%min non-zero
%min_vols = [curCGgrp.mVolCum];
%min_vols = cellfun(@(x) min(x), {curCGgrp.mVolCum});

% protect for unrealistically small min vols
min_vols = cellfun(@(x) min(x(x>0.1)), {curCGgrp.mVolCum});
vmin = min(min_vols);
vmin = round(vmin/vol_bin_step)*vol_bin_step;%round to nearest bin step
vol_bins = [log10(vmin) : vol_bin_step : log10(vmax)+vol_bin_step];
vol_bins = [0, 10.^vol_bins];



% % will write atlas for each fraction bin and for combined
% grp_flgs = [curCGgrp.mFxNum];
% 
% %unique fractionations
% fxs = unique(grp_flgs);
% % get number of entries for each unique fx number
% fx_counts = histc(grp_flgs,fxs);
% 
% % add combined fractionation schemes
% fxs = [-1 fxs];
% fx_counts = [sum(fx_counts) fx_counts];


cur_fx_counts = histc([curCGgrp.mFxNum],fxs(2:end));
fx_counts = [sum(cur_fx_counts) cur_fx_counts];
    
%% loop over each fractionation scheme (last entry is combined schemes)
for i=1:length(fxs),
    cur_fx = fxs(i);
    
    % get number of entries for each unique fx number

    cur_fx_pts = fx_counts(i);
    fx_flg = true(length(curCGgrp),1)';
    if ~(cur_fx<0), % cur_fx = -1; combines all fxs
        fx_flg = logical([curCGgrp.mFxNum]==cur_fx);
    end
    FxGrp = curCGgrp(fx_flg);
    FxComp = flgcomp(fx_flg);
    
    %% atlas generation
    
    AtlasTotal = zeros( length(dose_bins), length(vol_bins));
    AtlasComp = zeros( length(dose_bins), length(vol_bins));
    
    %% Rebin patient DVHs given dose_bins
    
    
    for ii=1:length(dose_bins),
        
        cur_dose = dose_bins(ii);
        cur_vols = zeros(cur_fx_pts,1); % patient volumes for current dose bin
        
        %% Loop over patients for given structure-outcome-fx
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
    
    FxAtlas(end,3:2:end) = {'nc(D(i),V(j),t(k))'};
    FxAtlas(end,4:2:end) = {'np(D(i),V(j),t(k))'};
    
       
    if cur_fx<0
       write_atlas_file_name = [atlas_file_name,'allfr']; 
    else
        write_atlas_file_name = [atlas_file_name,num2str(cur_fx),'fr'];
    end
    
   %atlas_file_name = [atlas_file_name,'_aoc.xls'];
    
    warning('off','MATLAB:xlswrite:AddSheet');
    disp(['Writing ',write_atlas_file_name,', sheet: ',cur_atlas_sheet_name,'...']);
    xlswrite(write_atlas_file_name,FxAtlas,cur_atlas_sheet_name,'A1');
    warning('on','MATLAB:xlswrite:AddSheet');
    

end
 %atlas location and file name
   
%     if fxs(i)<0, 
%         atlas_sheet_name = [atlas_sheet_name,'all'];
%     else
%         atlas_sheet_name = [atlas_sheet_name,num2str(fxs(i))];
%     end
%     

end


toc;
end