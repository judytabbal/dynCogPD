function results = go_decomposeConnectome_SS_ICA(state,varagin)

% This function applies ICA decomposition (JADE algorithm) to extract
% independent components using the concatenation approach of all temporal
% windows/trials/subjects to generate spatial maps and temporal signals of
% NCs independent components
%
% Input: 
%   state : structure containing the following field:
%         state.NCs: number of Results Independent Components (ICs)
%         state.conn: cell of length nb of subject, containing the path of cmat structure (output of compute_source_connectivity.m)
%         state.n_parcels: number of parcels or ROIs used
%
% Output:
%   results : structure containing results of Source Separation ICA Decomposition (connectivities maps and temporal signals)
%         results.NCs: number of Independent Components (ICs)
%         results.n_trials: number of total trials for the group
%         results.sub_trials: array of length (nb_sub) containing number of trials per subject
%         results.time: array of length (nb_temporalwindow) containing the time instant (in sec) of the center of each sliding window
%         results.signals: 3D matrix of dim [NCs,nb_temporalwindow,n_trials] 
%         results.maps: 3D matrix of dim [nROIs,nROIs,NCs] 


% 1. Concatenate all the data (temporal windows of all trials and subjects). 
n_subs      = length(state.conn);
cmat_all    = [];

disp('Loading and concatenating connectomes :')
ft_progress('init', 'text', 'Please wait...')
pause(0.01)
% loop across subjects
for ii = 1:n_subs
    
    cmat     = [];
    cmat_sub = [];
    cmat_2d  = [];
    cmat_red = [];   
    
    ft_progress(ii/n_subs, 'Processing connectome %d of %d', ii, n_subs);  
    load(state.conn{ii});
    
    % loop across trials
    for jj = 1:cmat.n_trials
        tmp = cmat.connectivity{jj};
        cmat_sub = cat(3,cmat_sub,tmp);
    end  
    
    % Do a DC correction on each connection
    cmat_sub = cmat_sub-repmat(mean(cmat_sub,3),[1 1 size(cmat_sub,3)]);   
    % flatten tensor into 2.5d (as the machine learning kids call it)
    cmat_2d = reshape(cmat_sub,state.n_parcels*state.n_parcels,length(cmat.time)*cmat.n_trials);
    % as each slice of the tensor is symmetric, remove redundant data to save some memory
    index = find(tril(squeeze(cmat_sub(:,:,1)))~=0);
    cmat_red = cmat_2d(index,:);     
    % concatenate onto the group level connectome data
    cmat_all = cat(2,cmat_all,cmat_red);   
    sub_trials(ii) = cmat.n_trials;   
    
end
disp('DONE')

% 2. Apply JADE-ICA algo
[B] = jader(cmat_all, state.NCs);
A=B'; 
sig=B*cmat_all; 

% 3. Post-process ICA results 
tmp = reshape(sig,state.NCs,size(sig,2)/sum(sub_trials),sum(sub_trials));
results.signals = tmp;
for ii = 1:state.NCs
    tmp = zeros(state.n_parcels,state.n_parcels);
    tmp(index) = A(:,ii);
    tmp = tmp + transpose(tmp);
    results.maps(:,:,ii) = tmp;
end

% 4. Complete the structure of results output
results.maps       = abs(results.maps);
results.NCs        = state.NCs;
results.n_trials   = sum(sub_trials);
results.sub_trials = sub_trials;
results.time       = round(cmat.time*100)/100;
