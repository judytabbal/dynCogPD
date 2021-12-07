 function [filters_grp,cmat_grp] = compute_source_connectivity(source,dFC,database,inparam,dataStruct)

% This function computes wMNE sources weights (=filters_grp) and PLV
% connectivity structure (=cmat_grp) of all subjects
%
% Inputs:
%   source: input structure of source
%   dFC: input structure of dynamic functional connectivity
%   database: input structure of database
%   inparam: input structure of useful input parameters
%   dataStruct : input structure of FIeldTrip data (output of create_dataStruct function)
%
% Outputs:
%   filters_grp: computed wMNE filter for the group (each cell represents one subj)
%   cmat_grp: path of computed dFC for the group (each cell represents one subj)


% Loop over all subjects
for i_sub = 1:database.nb_sub
    
    filters=[];
    cmat=[];
    
    % Source Reconstruction
    cfg             = [];
    cfg.prestim     = database.pre_samples/database.fs; % prestimulus in seconds for noise cov computation
    cfg.weightExp   = source.weightExp; % param for depth weighting
    cfg.weightLimit = source.weightLimit; % param for depth weighting limit
    cfg.SNR         = source.SNR; % Signal to Noise Ratio
    filters         = go_source_reconstruction(cfg,dataStruct{i_sub},inparam.subvol,inparam.subgrid,inparam.scout_scs.orients);
    
    filters_grp{i_sub} = filters;

    % Dynamic FC Computation
    cfg             = [];
    cfg.window.size = dFC.window.size; % sliding window length in seconds 
    cfg.window.step = dFC.window.step; % step between windows in seconds 
    cfg.bpfreq      = dFC.bpfreq; % frequency band of interest
    cfg.prestim     = database.pre_samples/database.fs; % prestimulus time in seconds
    cfg.poststim    = database.post_samples/database.fs; % poststimulus time in seconds
    cfg.conn_method = dFC.conn_method; % 'plv_dyn' for windowedPLV (only this option is implemented in code)
    cfg.labels      = inparam.scout_labels;
    cmat            = go_dynamicFC(cfg,filters,dataStruct{i_sub});
    
    cmat_grp{i_sub} = [inparam.path_conn '\cmat_PLV_' int2str(i_sub) '.mat']; 

    % Save cmat results
    save([inparam.path_conn '\cmat_PLV_' int2str(i_sub) '.mat'],'cmat');

end

