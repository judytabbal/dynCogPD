%% Make sure to open brainstorm software when running the following codes

%% Basic Params to Manipulate appropriately: RUN THIS SECTION ONCE IN THE BEGINNING

path_base = 'E:\DynamicEEG-Park'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mri_template = 'icbm'; %'icbm' for ICBM152
atlas = 'destrieux'; %'destrieux' for Destrieux148 atlas (same format as BrainStorm)
ROIs=148; %148 for destrieux
ind_atlas=2; %2 for destrieux
scs_mni=1; %1 for scs, 2 for mni coordinate system

% Add defaults fieldtrip functions (Run this section in the beginning)
tmp = which('ft_defaults');
if isempty(tmp)
    addpath([path_base '\Toolboxes\fieldtrip-20190224']); %add defaults fieldtrip
    ft_defaults
end


%% 1- Read MRI download from freesurfer  + Realign MRI 

mri=ft_read_mri([path_base '\Inputs\' mri_template '\icbm_avg_152_t1_tal_lin.nii']); 

cfg=[];
mri_realign=ft_volumerealign(cfg,mri);


%% 2- Scout Structure: destrieux with mri icbm (scs or mni coordinate systems)

load([path_base '\Inputs\' mri_template '\cortex_' mri_template '.mat']); %load the cortex 15002 template exported from Brainstorm 
load([path_base '\Inputs\' mri_template '\mri_' mri_template '_BS.mat']); %load the mri template exported from Brainstorm 

for i=1:ROIs
    % Get the vertex indices for the first scout of the atlas
    iSeed = cortex.Atlas(1,ind_atlas).Scouts(1,i).Seed;
    % Get the SCS coordinates for the corresponding vertices
    Vscs = cortex.Vertices(iSeed,:);
    if(scs_mni==1)
       centroids(i,:) = Vscs;
    elseif(scs_mni==2)
       centroids(i,:) = cs_convert(mri_BS, 'scs', 'mni', Vscs); %brainstorm should be opened to compute cs_convert
    end  
    iSeed = cortex.Atlas(1,ind_atlas).Scouts(1,i).Seed;
    Oscs = cortex.VertNormals(iSeed,:);
    orients(i,:)=Oscs;
end

sizeVertices=size(cortex.Vertices,1);
for j=1:sizeVertices
    allVscs = cortex.Vertices(j,:);
    if(scs_mni==1)
        allVert(j,:) = allVscs;
    elseif(scs_mni==2)
        allVert(j,:) = cs_convert(mri_BS, 'scs', 'mni', allVscs);
    end  
end

for i=1:ROIs  
    label{i,1}= cortex.Atlas(1,ind_atlas).Scouts(1,i).Label;    
end

scout_labels=label;
if(scs_mni==1)
    scout_scs.centroids=centroids;
    scout_scs.faces=cortex.Faces;
    scout_scs.vertices=allVert;
    scout_scs.orients=orients;
elseif(scs_mni==2)
    scout_mni.centroids=centroids;
    scout_mni.faces=cortex.Faces;
    scout_mni.vertices=allVert;
    scout_mni.orients=orients;
end


%% 3- Source Grid (fieldtrip format)

%Load Atlas SCS (CTF)
load([path_base '\Inputs\' mri_template '\' atlas '\' atlas '_scs_' mri_template '.mat']);

cfg                 = [];
cfg.grid.pos        = scout_scs.centroids*1000; %scs/mm
cfg.grid.inside     = ones(size(scout_scs.centroids,1),1);
subgrid             = ft_prepare_sourcemodel(cfg);


%% 4- Electrode format conversion from Brainstorm to fieldTrip

%open Brainstorm or insert to the path
load([path_base '\Inputs\' mri_template '\elec_bst.mat']); %electrode structure extracted from BS
[ftElec, ftGrad] = out_fieldtrip_channel(elec_bst);
elec_BS_mm=ft_convert_units(ftElec,'mm');

% % visualize electrodes alignment (elec_BS_mm) with headmodel (subvol)
% figure();
% ft_plot_sens(elec_BS_mm,'facecolor','r','elecsize',10); 
% hold on
% ft_plot_vol(subvol, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;


%% 5- Example of indata variable construction based on your own data 

% 'indata.mat' represents the set of all preprocessed EEG trials for all
% subjects. This data should be loaded in the 'run_pipeline.m' code for
% further segmentation and to appropriatly create data strcuture in FieldTrip.

% Note: the variable parameters below are not related to the real data used,
% they are chosen randomly to illustrate how should be the 'indata.mat' format

nb_sub=10; % number of subjects
nb_trials_persub=[2 4 5 3 8 5 4 9 10 6]; % number of trials for each subject
nb_channels=32; % number of channels (electrodes)
nb_samples=125; % number of total samples in trial

for i=1:nb_sub 
    for j=1:nb_trials_persub(i) 
        indata{i}{j}=randn(nb_channels,nb_samples); % randn is used only as example, however you should replace it with your own data
    end
end
