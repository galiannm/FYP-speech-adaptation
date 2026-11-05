
%% ============SETUP============
clear all
close all

addpath ../../libs/cnsp_utils
addpath ../../libs/cnsp_utils/cnd
addpath ../../libs/mTRF-Toolbox_v2/mtrf
addpath ../../libs/NoiseTools
addpath ../../libs/eeglab
eeglab

%% ============SELECTING DATASETS & STIMULUS FEATURES============

% ----Envelope----
% dataMainFolders = {'AAD_KULeuven/', 'AliceSpeech/', 'diliBach/', 'FDSpeech_L1/', 'FDSpeech_L2/', ...
%     'LalorNatSpeech/', 'LalorNatSpeechReverse/', 'SparrKULee1', 'childStoriesEEG_Sysoeva/', 'cocktailAttSwitch/', ...
%     'cocktailAttSwitch/', 'fNIRS_podcast/', 'musicImagery/', 'polyphonicBach/', 'spectrotemporalModulationsIgnacio/', ...
%     'hyperscanningAmir_Listening/', 'podcastListeningEEG/', 'ComedyShowEEG/', 'speechEmotionAdultEEG/', 'speechTrustEEG/', ...
%     'sungSpeechCantisani_melody', 'sungSpeechCantisani_song'};

%     stimFeatureDataset = [1, 1, 1, 1, 1, ...
%                             1, 1, 1, 1, 1, ...
%                             2, 1, 1, 3, 1, ...
%                             1, 1, 1, 1, 1, ...
%                             1, 1];

% featureName = "envelope"

% ----Word Onset----
dataMainFolders = {'AliceSpeech/', 'FDSpeech_L1/', 'LalorNatSpeech/', 'childStoriesEEG_Sysoeva/', ...
    'cocktailAttSwitch/', 'cocktailAttSwitch/', 'fNIRS_podcast/', 'hyperscanningAmir_Listening/', ...
    'podcastListeningEEG/', 'ComedyShowEEG/', 'speechEmotionAdultEEG/', 'speechTrustEEG/', };

stimFeatureDataset = [2, 5, 4, 2, 3, 4, 3, 4, 3, 11, 3, 3];

featureName = "wordOnset";

toPreprocess = zeros(1, length(dataMainFolders));

%% Dataset Paths & Output
scriptFolder = fileparts(mfilename('fullpath'));
dataCNDSubFolder = 'dataCND/';

for iiDataset = 1:length(dataMainFolders)
    datasetName = dataMainFolders{iiDataset};
    fprintf('\nProcessing dataset: %s\n', datasetName);

    dataMainFolder = fullfile(scriptFolder, '..', '..', 'datasets2', dataMainFolders{iiDataset});
    dataCNDfolder = fullfile(dataMainFolder, 'dataCND');

    % with 'Avg' we compute the mean voltage across the electrodes
    reRefType = 'Avg'; % or 'Mastoids'
    bandpassFilterRange = [1,8]; % Hz

    downFs = [];
    eegFolder = fullfile(dataCNDfolder, '1-8Hz');
    fprintf('Looking in: %s\n', eegFolder);
    fprintf('Dataset folder exists: %d\n', isfolder(dataMainFolder));
    fprintf('dataCND folder exists: %d\n', isfolder(dataCNDfolder));
    fprintf('1-8Hz folder exists: %d\n', isfolder(eegFolder));
    eegFilenames = dir(fullfile(eegFolder, 'pre_dataSub*.mat')); % participants files
    nSubs = length(eegFilenames);
    fprintf('Found %d subjects\n', nSubs);

    %% LOAD AND INSPECT EEG
    % Check that EEG files were found
    if nSubs == 0
        warning('No EEG files found for %s', datasetName);
        continue
    end
    
    % Load the first subject as an example
    eegFilename = fullfile(eegFolder, eegFilenames(1).name);
    
    fprintf('Loading EEG file: %s\n', eegFilenames(1).name);
    
    load(eegFilename, 'eeg');
    
    % Display the fields stored in the EEG structure
    disp('EEG structure fields:');
    disp(fieldnames(eeg));
    
    % Display sampling frequency
    fprintf('Sampling frequency: %g Hz\n', eeg.fs);
    
    % Display number of trials
    nTrials = numel(eeg.data);
    fprintf('Number of trials: %d\n', nTrials);
    
    % Display dimensions of the first trial
    [nSamples, nChannels] = size(eeg.data{1});
    
    fprintf('First trial: %d samples, %d channels\n', ...
        nSamples, nChannels);
    
    % Calculate the duration of the first trial
    duration = nSamples / eeg.fs;
    
    fprintf('First trial duration: %.2f seconds\n', duration);

    %% ============INSPECT STIMULUS FILES============

    %% Load Stim Info
    % Check whether there is one stimulus file for all subjects
    sharedStimFile = fullfile(dataCNDfolder, 'dataStim.mat');

    % All subjects share the same stimulus file (dataStim.mat)
    if isfile(sharedStimFile)
        singleDataStim = true;
        stimFilename = sharedStimFile;
    % All subjects have a seperate stimulus file (dataStim*.mat)
    else
        singleDataStim = false;
        stimFiles = dir(fullfile(dataCNDfolder, 'dataStim*.mat'));

        if isempty(stimFiles)
            error('No stimulus files found for %s', datasetName);
        end

        stimFilename = fullfile(dataCNDfolder, stimFiles(1).name);

    end
    load(stimFilename, 'stim');
    fprintf('Loaded stimulus file: %s\n', stimFilename);

    %% INSPECT STIMULUS STRUCTURE

    disp('Stimulus fields:');
    disp(fieldnames(stim));
    
    disp('Available stimulus features:');
    disp(stim.names);
    
    fprintf('Stimulus sampling frequency: %g Hz\n', double(stim.fs));
    fprintf('Number of stimulus trials: %d\n', numel(stim.data));
    
    disp('Dimensions of first stimulus trial:');
    disp(size(stim.data{1}));

    %% INSPECT EXPERIMENTAL CONDITIONS
    
    % Number of stimulus trials
    nStimTrials = size(stim.data, 2);
    
    fprintf('Number of stimulus trials: %d\n', nStimTrials);
    
    % Check whether condition indices exist
    if isfield(stim, 'condIdxs') && ~isempty(stim.condIdxs)
    
        disp('Experimental condition indices:');
        disp(stim.condIdxs);
    
        disp('Unique conditions:');

        if iscell(stim.condIdxs)
            conditionValues = cell2mat(stim.condIdxs);
        else
            conditionValues = stim.condIdxs;
        end
        
        disp(unique(conditionValues, 'stable'));
    
    else
    
        disp('No condition indices found.');
        disp('Assuming all trials belong to condition 1.');
    
        stim.condIdxs = ones(1, nStimTrials);
    
    end
end