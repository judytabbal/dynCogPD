function dataStruct = create_dataStruct(database,inparam)

% This function calculates data structure (fieldTrip format) from input database parameters as follows:
%
% Inputs:
%   database : structure with the following fields:
%              database.nb_sub       : number of subjects
%              database.nb_trials    : array of length (nb_sub) that contains number of trials per subject 
%              database.nb_chan      : number of channels (electrodes)
%              database.fs           : sampling frequency
%              database.onset        : sample index for onset 
%              database.pre_samples  : nb samples before trial onset 
%              database.post_samples : nb samples after trial onset 
%              database.indata       : preprocessed EEG trials for all subjects, as a cell of length (nb_sub), each cell is 2D matrix of dim: [nb_chan,nb_samples]
%   inparam  : strcuture, should include elec_BS_m field
%
% Output: cell of length nb_sub, each cell is a data strcuture of fieldtrip format
%   dataStruct.fsample: int indicating sampling frequency of data (in our case = fs)
%   dataStruct.elec: fieldtrip structure of electrodes (in our case = elec_BS_mm)
%   dataStruct.trial: the segmented data, cell of length nb_trials, each cell is 2D matrix [nb_chan*nb_samples]
%   dataStruct.time: the time in sec from -prestim to +poststim, cell of length nb_trials, each cell is 1D matrix [1*nb_samples]
%   dataStruct.label: cell of length nb_chan, each cell is a string of channels or electrodes labels (for example in case of EGI257, we can label electrodes from E1 to E256 and Cz for the last electrodes as a reference)
%
% DATA TRIALS SHOULD BE PREPROCESSED


% Loop over subjects
for i_sub=1:database.nb_sub

    % Build dataStruct fields for each subject
    data         = [];
    data.fsample = database.fs; 
    data.elec    = inparam.elec_BS_mm;

    for j=1:database.nb_trials(i_sub)
        % Segment trial relative to onset 
        Bepoched=database.indata{i_sub}{j}(:,database.onset-database.pre_samples:database.onset+database.post_samples);
        
        % Centralize (to zero) segmented trial
        Bepoched_demean=Bepoched-mean(Bepoched,2); 
        
        data.trial{j}= Bepoched_demean;
        data.time{j}=[-database.pre_samples:1:database.post_samples]/database.fs; 
    end
 
    label={};
    for i_chan=1:database.nb_chan-1
        label{i_chan}=['E' int2str(i_chan)];
    end
    label{end+1}='Cz';  
    data.label=label;

    dataStruct{i_sub}=data;
end

