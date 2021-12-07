function filters=go_source_reconstruction(opt,data,subvol,subgrid,sources_orients)


% 1. Calculate Lead fields (Forward Model)
cfg             = [];
cfg.headmodel   = subvol;
cfg.elec        = data.elec;
cfg.grid.pos    = subgrid.pos;
cfg.normalise   = 'yes';
cfg.rankreduce  = 3; 
lf              = ft_prepare_leadfield(cfg); 

   
% 2. Timelock for noise covariance estimation
cfg                     = [];
cfg.covariance          = 'yes';
cfg.window              = [-opt.prestim 0]; %baseline for noise cov (computed from prestim in seconds to 0)
tlk_noise                = ft_timelockanalysis(cfg,data);
noise_cov=tlk_noise.cov;

% 3. Compute wMNE filters
filters = ComputeWMNE(noise_cov,cell2mat(lf.leadfield),lf.pos,sources_orients,opt.weightExp,opt.weightLimit,opt.SNR);
