function runGenerateDvhActuarialAtlas

analy_str='cwp';

%% TODO Pooled analysis
% need to incorporate institutes, somehow
if isequal(analy_str,'pld'),

    %structure_strs={'lungs'};
    %toxicity_strs={'rp'};
    %data_loc_str='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta_comb.mat'
    %atlas_loc_str='Z:/elw/MATLAB/meta_analy/atlases';
    %GenerateDvhAtlas('pld','lungs','rp','

elseif isequal(analy_str,'nfz'),
    
    %% NFZ Analysis
    structure_strs = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
    %toxicity_strs = {'rp','pultox','esotox'};
    toxicity_strs = {'pultox'};
    %structure_strs = {'ILUNG'};
    %toxicity_strs = {'pultox'};
    
    data_loc_str='Z:/elw/MATLAB/nfz_analy/meta_data/';
    atlas_loc_str='Z:/elw/MATLAB/nfz_analy/atlases/latest/';

elseif isequal(analy_str,'bpx');
    %% Bpx Analysis
    structure_strs = {'bp'};
    toxicity_strs = {''};
    
    data_loc_str='Z:/elw/MATLAB/bpx_analy/meta_data/BPx_DiVj_DVHs_fx-1_a2bInf.mat';
    atlas_loc_str='Z:/elw/MATLAB/bpx_analy/atlases/';
    
elseif isequal(analy_str,'cwp'),
    
    %% NFZ Analysis
    structure_strs = {'cw'};
    toxicity_strs = {''};
    
    data_loc_str='Z:/elw/MATLAB/cw_analy/meta_data/MUTTER_MASTER_ChestWall_Cox_DiVj_DVHs_fx-1_a2bInf.mat';
    atlas_loc_str='Z:/elw/MATLAB/cw_analy/atlases/latest/';
    
else
    error(['runGenerateDvhAtlas:: Unknown analy_str: ',analy_str]);
end

%% Run
for i=1:length(structure_strs)
    cur_structure_str = structure_strs{i};
    for j=1:length(toxicity_strs);
        cur_toxicity_str = toxicity_strs{j};
        
        GenerateDvhActuarialAtlas(analy_str,cur_structure_str,cur_toxicity_str,...
                         data_loc_str,atlas_loc_str);
    end
end
end