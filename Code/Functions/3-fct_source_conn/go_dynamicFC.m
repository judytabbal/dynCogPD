function cmat=go_dynamicFC(opt,filters,data)

% Inputs:
%   filters : source weights 
%   data    : fieldTrip data structure of one subject
%   opt     : structure containing the following fields:
%             opt.window.size : sliding window length in sec
%             opt.window.step : step between windows in sec
%             opt.bpfreq      : frequency band of interest, ex [12 25]
%             opt.prestim     : prestimulus time in sec
%             opt.poststim    : poststimulus time in sec
%             opt.conn_method : 'plv_dyn' for windowedPLV (only this option is implemented in code)
%             opt.labels      : scout labels
%
% Outputs:
%   cmat: dFC structure of one subject containing the following fields:
%         cmat.connectivity : cell of length nb of trials, each is 3D matrix of dim [nROIs,nROIs,nb_windows]
%         cmat.time : array of length nb_windows, containing the time of center of each temporal window
%         cmat.n_trials : number of trials
%         cmat.window : window specification (size, step)
%         cmat.bpfreq : ex [12 25]
%         cmat.n_parcels : number of ROIs
%         cmat.conn_type: connectivity type: 'plv_dyn'



% 1. Define general variables
n_parcels = size(filters,1);
opt.trial_length=abs(opt.prestim)+abs(opt.poststim);
tmp = (0+opt.window.size/2):opt.window.step:(opt.trial_length-opt.window.size/2);
n_windows = length(tmp); % number of windows per trial
n_trials  = length(data.trial);
n_samples = floor(opt.window.size * data.fsample);
n_shifts  = floor(opt.window.step * data.fsample);
ti = 1+(0:(n_windows-1))*n_shifts; % start sample index of each window
tf = n_samples+(0:(n_windows-1))*n_shifts; % end sample index of each window

% 2. Build the cmat structure
cmat                    = [];
disp('Generating connectivity matrices:')
ft_progress('init', 'text', 'Please wait...')

% loop through all trials to generate connectomes
for ii = 1:n_trials
    ft_progress(ii/n_trials, 'Processing trial %d from %d', ii, n_trials);
    VE = [];
    VE.raw = filters*data.trial{ii};
    data_input{ii}=VE.raw;
end

eegData=zeros(n_parcels,size(data_input{1},2),n_trials);
for tr=1:n_trials
    eegData(:,:,tr)=data_input{tr};
end

fmin=opt.bpfreq(1);
fmax=opt.bpfreq(2);

tmp2 = data.time{1};
for ii = 1:n_windows
    time_wind(ii) = mean(tmp2(ti(ii):tf(ii)));
end

ROIs_labels=opt.labels;

% 3. Calculate dynamic PLV 
[plvs,N]=go_CalculateDynWindPLV(data_input,data.fsample,fmin,fmax,opt.window.size,opt.window.size-opt.window.step,opt.window.step,n_windows);

for t=1:n_trials
    plvs_perm{t}=permute(plvs{t},[2 3 1]);
    for w=1:n_windows
        plvs_sym{t}(:,:,w)=squeeze(plvs_perm{t}(:,:,w))+squeeze(plvs_perm{t}(:,:,w))';
    end
    cmat.connectivity{t}=plvs_sym{t};
end


ft_progress('close')
disp('DONE')

% 4. Complete cmat structure infos
cmat.time      = time_wind;
cmat.n_trials  = n_trials;
cmat.window    = opt.window;
cmat.bpfreq    = opt.bpfreq;
cmat.n_parcels = n_parcels;
cmat.conn_type = opt.conn_method;
