
%% ============SETUP============
clear all
close all

addpath ../libs/cnsp_utils
addpath ../libs/cnsp_utils/cnd
addpath ../libs/mTRF-Toolbox_v2/mtrf
addpath ../libs/NoiseTools
addpath ../libs/eeglab
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

%% Dataset Paths & Output
scriptFolder = fileparts(mfilename('fullpath'));
dataCNDSubFolder = 'dataCND/';
resultDir = fullfile(scriptFolder, 'resultsChunks', featureName);
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

for iiDataset = 1:length(dataMainFolders)
    datasetName = dataMainFolders{iiDataset};
    fprintf('\nProcessing dataset: %s\n', datasetName);

    dataMainFolder = fullfile(scriptFolder, '..', 'datasets2', dataMainFolders{iiDataset});
    dataCNDfolder = fullfile(dataMainFolder, 'dataCND');

    % with 'Avg' we compute the mean voltage across the electrodes
    reRefType = 'Avg'; % or 'Mastoids'
    bandpassFilterRange = [1,8]; % Hz

    downFs = [];
    eegFolder = fullfile(dataCNDfolder, '1-8Hz');
    eegFilenames = dir(fullfile(eegFolder, 'pre_dataSub*.mat')); % participants files
    nSubs = length(eegFilenames);
    fprintf('Found %d subjects\n', nSubs);

    %% ============The multivariate Temporal Response Function============
    % 0-10s, 10-20s, etc

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

    % Downsampling stim if necessary
    stim.data = cellfun(@double, stim.data, 'UniformOutput', false);
    stim.fs = double(stim.fs);
    downFs = double(downFs);

    if ~isfield(stim,'condIdxs')
        stim.condIdxs = ones(size(stim.data,2),1);
    end

    %% TRF Parameters and time windows
    % hyperparameters
    tmin = -100;
    tmax = 600; % 400 for envelope
    lambdas = [1e-6,1e-3,1e-4,1e-2,1e-1,1e0,1e2,1e3,1e4]; % regularisation param
    dirTRF = 1; % Forward TRF model (direction)

    stimIdx = stimFeatureDataset(iiDataset);  % stim feature for the current dataset

    chunkLength = 10
    maxChunks = 20

    for cond2select = 1:length(unique(stim.condIdxs))
        for ii = 1:maxChunks
            fprintf('Condition %d, chunk %d\n', cond2select, ii);

            % TRF
            rAll = []; % one prediction correlation per subject (averaged across electrodes)
            rAllElec = []; % prediction correlation for each electrode and subject
            modelAll = struct([]); % TRF model fitted for each subject
            includedSubjectIdx = []; % Ids for subject successfully included
            chanlocs = [];

            wordCountMeanAll = [];
            wordCountStdAll = [];   
            
            %% Loading each subject's EEG and stimulus data
            for sub = 1:nSubs
                toBreak = false;

                subIdx = regexp(eegFilenames(sub).name, '\d+', 'match');
                subIdx = str2double(subIdx{1});

                fprintf('\nProcessing subject %d\n', subIdx);

                eegFilename = fullfile(eegFolder, eegFilenames(sub).name);
                load(eegFilename, 'eeg');

                if ~singleDataStim
                    stimFilename = fullfile(dataCNDfolder, sprintf('dataStim%d.mat', subIdx));
                    load(stimFilename, 'stim');
                    stim.fs = double(stim.fs);
                end

                stimFeature = stim;

                % select feature name
                if iscell(stimFeature.names)
                    stimFeature.names = stimFeature.names{stimIdx};
                end
                stimFeature.data = stimFeature.data(stimIdx, :);
                if ~isfield(stimFeature, 'condIdxs') || isempty(stimFeature.condIdxs)
                    stimFeature.condIdxs = ones(1, numel(stimFeature.data));
                end

                %% Count words per interval
                wordCounts = 0;
                if strcmpi(featureName, 'wordOnset')
                    originalFs = double(stimFeature.fs);
                    wordStart = 1 + round(originalFs * chunkLength * (ii - 1));
                    wordEnd = round(originalFs * chunkLength * ii);
                    wordCounts = nan(1, numel(stimFeature.data));

                    for tr = 1:numel(stimFeature.data)
                        wordSignal = stimFeature.data{tr};
                        if wordEnd <= size(wordSignal, 1)
                            wordChunk = wordSignal(wordStart:wordEnd, :);
                            wordCounts(tr) = nnz(wordChunk);
                        end
                    end
                end

                % Match stimulus sampling freq to preprocessed EEG
                eeg.fs = double(eeg.fs);
                stimFeature.fs = double(stimFeature.fs);

                if stimFeature.fs ~= eeg.fs

                    fprintf('Resampling STIM from %g Hz to %g Hz\n', stimFeature.fs, eeg.fs);

                    oldFs = stimFeature.fs;
                    newFs = eeg.fs;

                    for tr = 1:numel(stimFeature.data)
                        stimFeature.data{tr} = resample(double(stimFeature.data{tr}), newFs, oldFs);
                    end

                    stimFeature.fs = newFs;
                end

                assert(stimFeature.fs == eeg.fs, 'EEG and STIM sampling frequencies do not match');

                % Save actual sampling frequency
                downFs = eeg.fs;
                fprintf('EEG fs = %g Hz | STIM fs = %g Hz\n', eeg.fs, stimFeature.fs);

                % Skip subjects where EEG and stimulus trial counts do not match
                if numel(stimFeature.data) ~= numel(eeg.data)
                    warning('Skipping subject %d: STIM has %d trials, EEG has %d trials', subIdx, numel(stimFeature.data), numel(eeg.data));
                    continue
                end

                [stimFeature,eeg] = cndCheckStimNeural(stimFeature,eeg);
                stimFeature = cndNormalise(stimFeature);
                eeg = cndNormalise(eeg);
                [stimFeature,eeg] = cndCutSameLength(stimFeature,eeg);

                %% TRF crossvalidation
                sampStart = 1+eeg.fs*chunkLength*(ii-1);
                sampEnd = eeg.fs*chunkLength*(ii);

                for jjTr = 1:length(eeg.data)
                    if sampEnd > size(eeg.data{jjTr}, 1)
                        toBreak = true;
                    end
                end
                if toBreak
                    fprintf('Skipping subject %d: insufficient data for chunk %d\n', subIdx, ii);
                    continue
                end
                % extracting the data
                stimFeature.data = cellfun(@(x) x(sampStart:sampEnd,:), stimFeature.data,'UniformOutput',false);
                eeg.data = cellfun(@(x) x(sampStart:sampEnd,:), eeg.data,'UniformOutput',false);
                
                % Cross validation & choosing the best lambda
                disp('Running mTRFcrossval')

                conditions = unique(stimFeature.condIdxs, 'stable');

                % Skip subjects without the requested condition
                if cond2select > numel(conditions)
                    warning('Skipping subject %d: condition %d unavailable', subIdx, cond2select);
                    continue
                end

                if iscell(stimFeature.condIdxs)
                    trials2keep = strcmp(stimFeature.condIdxs, conditions{cond2select});
                else
                    trials2keep = stimFeature.condIdxs == conditions(cond2select);
                end

                % cross validation
                [stats,t] = mTRFcrossval( ...
                    stimFeature.data(trials2keep), ... % selected trails
                    eeg.data(trials2keep), ... % corresponding EEG trials
                    eeg.fs, ...    % sampling freq
                    dirTRF, ...    % forward model
                    tmin,tmax, ... % time lag
                    lambdas, ...   % regularization values to evaluate
                    'verbose',0);

                % Calculating optimal lambda
                [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));
                disp(['r = ',num2str(maxR)])

                rAll(end+1) = maxR;
                rAllElec(:,end+1) = squeeze(mean(stats.r(:,bestLambda,:),1));

                %% Fit TRF model ith optimal regularisation parameter
                disp('Running mTRFtrain')

                conditions = unique(stimFeature.condIdxs, 'stable');

                if iscell(stimFeature.condIdxs)
                    trials2keep = strcmp(stimFeature.condIdxs, conditions{cond2select});
                else
                    trials2keep = stimFeature.condIdxs == conditions(cond2select);
                end
                
                % fitting the final TRF model
                model = mTRFtrain( ...
                    stimFeature.data(trials2keep), ... % signals from selected trials
                    eeg.data(trials2keep), ... % eeg signals from those trials
                    eeg.fs, ... % sampling freq
                    dirTRF, ... % forward model
                    tmin, tmax, ... % time lag
                    lambdas(bestLambda), ... % regularisation
                    'verbose', 0); % disale detailed output
                % output: 
                % - model.w = estimated TRF weights
                % - model.t = time lags
                % - model.b = intercept

                %% Store TRF model
                if isempty(modelAll)
                    modelAll = model;

                    if isfield(eeg,'chanlocs') && ~isempty(eeg.chanlocs)
                        chanlocs = eeg.chanlocs;
                    else
                        chanlocs = [];
                    end
                else
                    modelAll(end+1) = model;
                end

                includedSubjectIdx(end+1) = subIdx;
                disp(['Mean r = ',num2str(mean(rAll))])
                
                selectedWordCounts = wordCounts(trials2keep);
                meanWords = mean(selectedWordCounts, 'omitnan');
                stdWords = std(selectedWordCounts, 0, 'omitnan');
                
                wordCountMeanAll(end+1) = meanWords;
                wordCountStdAll(end+1) = stdWords;

                fprintf('Mean words per trial: %.2f +/- %.2f\n', meanWords, stdWords);
            end
            
            %% Saving the TRF results
            if ~isempty(modelAll)
                filename = sprintf( ...
                    'trfOut_10sChunk_%s_cond%d_stim%d_%d.mat', ...
                    char(strip(string(dataMainFolders{iiDataset}), '/')), ...
                    cond2select, stimIdx, ii);

                save(fullfile(resultDir, filename), ...
                    'chanlocs', 'modelAll', 'rAllElec', ...
                    'bandpassFilterRange', 'downFs', 'reRefType', ...
                    'dataMainFolder', 'nSubs', 'stimIdx', ...
                    'tmin', 'tmax', 'lambdas', 'includedSubjectIdx', ...
                    'wordCountMeanAll', 'wordCountStdAll');
            end
        end
    end
end