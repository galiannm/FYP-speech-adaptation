% Cognition and Natural Sensory Processing (CNSP) Workshop
% Forward TRF - trial-by-trial envelope analysis
%
% Based on the CNSP example script.
% The original modelling logic is preserved:
%   - select one trial for mTRFcrossval (3 splits)
%   - select the best lambda for that trial
%   - train mTRFtrain on ALL trials with that lambda
%
% This version loads already preprocessed EEG and supports:
%   dataStim.mat OR participant-specific dataStimN.mat
%   origTrialPosition where present
%   different numbers of trials across participants

clear all
close all

% Run from the CNSP directory (the one containing datasets2 and libs).
addpath('libs/cnsp_utils')
addpath('libs/cnsp_utils/cnd')
addpath('libs/mTRF-Toolbox_v2/mtrf')
addpath('libs/NoiseTools')
addpath('libs/eeglab')
% eeglab  % Uncomment only if needed to initialise EEGLAB.

%% Datasets and envelope features (same ordering as requested)
dataMainFolders = { ...
    'AAD_KULeuven', 'AliceSpeech', 'diliBach', 'FDSpeech_L1', 'FDSpeech_L2', ...
    'LalorNatSpeech', 'LalorNatSpeechReverse', 'SparrKULee1', ...
    'childStoriesEEG_Sysoeva', 'cocktailAttSwitch', 'cocktailAttSwitch', ...
    'fNIRS_podcast', 'musicImagery', 'polyphonicBach', ...
    'spectrotemporalModulationsIgnacio'};

stimFeatureDataset = [1,1,1,1,1, 1,1,1,1,1, 2,1,1,3,1];

%% Settings
rootFolder = 'datasets2';  % Your current MATLAB folder is CNSP.
preprocessedBand = '1-8Hz'; % Use preprocessed EEG from this subfolder.
resultsFolder = 'resultsTrials_envelope';
makePlots = false;         % Enable after a successful run; saves one PNG per trial.

% IMPORTANT: data are already preprocessed. No preprocessing is repeated.
bandpassFilterRange = [1,8];
tmin = -100;
tmax = 400;
lambdas = [1e-6,1e-3,1e-4,1e-3,1e-2,1e-1,1e0,1,1e2,1e3,1e4];
dirTRF = 1; % Forward TRF model
% Be careful: backward models with many electrodes and large time-windows
% can require long computational times.

if ~isfolder(rootFolder)
    error('Cannot find %s. Set MATLAB Current Folder to the CNSP directory.',rootFolder);
end
if ~exist(resultsFolder,'dir'), mkdir(resultsFolder); end

%% Loop through dataset / feature combinations
for iiDataset = 1:numel(dataMainFolders) % TO CHANGE BACK TO 1
    datasetName = dataMainFolders{iiDataset};
    stimIdx = stimFeatureDataset(iiDataset);
    fprintf('\n*** %s (feature %d) ***\n',datasetName,stimIdx);

    dataCNDPath = fullfile(rootFolder,datasetName,'dataCND');
    eegFolder = fullfile(dataCNDPath,preprocessedBand);
    if ~isfolder(dataCNDPath)
        warning('Missing dataCND folder: %s. Skipping.',dataCNDPath);
        continue
    end

    % Load existing 1-8Hz preprocessed EEG; do not re-filter or re-reference.
    eegFilenames = dir(fullfile(eegFolder,'pre_dataSub*.mat'));
    if isempty(eegFilenames)
        warning('No preprocessed EEG found in %s. Skipping.',eegFolder);
        continue
    end
    nSubs = numel(eegFilenames);

    % Discover whether the dataset uses a shared or participant-specific stim.
    sharedStimFilename = fullfile(dataCNDPath,'dataStim.mat');
    hasIndividualStim = ~isfile(sharedStimFilename);
    if hasIndividualStim && isempty(dir(fullfile(dataCNDPath,'dataStim*.mat')))
        warning('No stimulus files found in %s. Skipping.',dataCNDPath);
        continue
    end

    % Establish the maximum trial count, including subject-specific lengths.
    subjectIds = nan(1,nSubs);
    subjectStimFiles = cell(1,nSubs);
    subjectNTrials = zeros(1,nSubs);
    for sub = 1:nSubs
        token = regexp(eegFilenames(sub).name,'pre_dataSub(\d+)\.mat$','tokens','once');
        if isempty(token)
            error('Cannot identify subject from filename: %s',eegFilenames(sub).name);
        end
        subjectIds(sub) = str2double(token{1});
        if hasIndividualStim
            subjectStimFiles{sub} = fullfile(dataCNDPath, ...
                sprintf('dataStim%d.mat',subjectIds(sub)));
        else
            subjectStimFiles{sub} = sharedStimFilename;
        end
        if ~isfile(subjectStimFiles{sub})
            warning('Missing stimulus for %s; subject will be skipped.',eegFilenames(sub).name);
            continue
        end
        S = load(subjectStimFiles{sub},'stim');
        subjectNTrials(sub) = size(S.stim.data,2);
    end
    nTr = max(subjectNTrials);
    if nTr == 0, continue; end

    outputDatasetFolder = fullfile(resultsFolder,sprintf('%s_stim%d',datasetName,stimIdx));
    if ~exist(outputDatasetFolder,'dir'), mkdir(outputDatasetFolder); end

    %% Loop through the stimulus trial numbers
    for ii = 1:nTr
        fprintf('\n%s | feature %d | stimulus trial %d/%d\n', ...
            datasetName,stimIdx,ii,nTr);

        % Keep results indexed by subject file, even if a trial is missing.
        rAll = nan(1,nSubs);
        bestLambdaAll = nan(1,nSubs);
        trialIdxAll = nan(1,nSubs);
        presentationPositionAll = nan(1,nSubs);
        rAllElec = []; % electrodes x subjects, filled with NaN for missing subjects
        modelAll = struct([]); % models for valid participants (same model type)
        modelSubjectIds = [];
        chanlocs = [];
        downFs = NaN;
        nValid = 0;
        plotEEG = [];

        for sub = 1:nSubs
            if subjectNTrials(sub) < ii || ~isfile(subjectStimFiles{sub})
                continue
            end

            %% Loading preprocessed EEG data
            eegPreFilename = fullfile(eegFolder,eegFilenames(sub).name);
            fprintf('Subject %d: loading %s\n',subjectIds(sub),eegFilenames(sub).name);
            E = load(eegPreFilename,'eeg');
            eeg = E.eeg;

            %% Loading Stimulus data
            S = load(subjectStimFiles{sub},'stim');
            stim = S.stim;
            % Downsampling stim if necessary (match EEG sampling rate).
            %% Resample stimulus to match EEG sampling rate
            if stim.fs ~= eeg.fs
            
                fprintf('Resampling stimulus from %g Hz to %g Hz\n', ...
                    stim.fs, eeg.fs);
            
                originalFs = double(stim.fs);
                targetFs = double(eeg.fs);
            
                [p,q] = rat(targetFs/originalFs,1e-10);
            
                for tr = 1:numel(stim.data)
            
                    % Samples x features
                    if isempty(stim.data{tr})
                        continue;
                    end
                
                    stim.data{tr} = double(stim.data{tr});
                    stim.data{tr} = resample(stim.data{tr},p,q);
            
                end
            
                stim.fs = targetFs;
            
            end
            if stim.fs ~= eeg.fs
                warning('%s subject %d: stim fs %g differs from EEG fs %g; skipping.', ...
                    datasetName,subjectIds(sub),stim.fs,eeg.fs);
                continue
            end
            if ~isfield(stim,'condIdxs')
                stim.condIdxs = ones(size(stim.data,2),1);
            end
            if stimIdx > size(stim.data,1)
                warning('%s subject %d: feature %d does not exist; skipping.', ...
                    datasetName,subjectIds(sub),stimIdx);
                continue
            end

            %% Selecting feature of interest ('stimIdx' feature)
            stimFeature = stim;
            if iscell(stimFeature.names)
                stimFeature.names = stimFeature.names{stimIdx};
            elseif isstring(stimFeature.names)
                stimFeature.names = stimFeature.names(stimIdx);
            elseif ischar(stimFeature.names)
                stimFeature.names = stimFeature.names;
            else
                warning('Unsupported feature-name type: %s', class(stimFeature.names));
                stimFeature.names = sprintf('Feature_%d',stimIdx);
            end
            stimFeature.data = stimFeature.data(stimIdx,:); % envelope

            %% Making sure that stim and neural data have the same length
            % The trial may end a few seconds after the end of the audio
            % e.g., the neural data may include the break between trials
            % It would be best to do this chunking at preprocessing, but let's
            % check here, just to be sure
            [stimFeature,eeg] = cndCheckStimNeural(stimFeature,eeg);

            %% Standardise stim data (preserving the ratio between features)
            % This is thought for continuous signals e.g., speech envelope, eeg
            stimFeature = cndNormalise(stimFeature);

            %% Standardise neural data (preserving the ratio between channels)
            eeg = cndNormalise(eeg);
            
            %% Check that stimulus and EEG trials correspond

            nStimTrials = size(stimFeature.data,2);
            nEEGTrials = size(eeg.data,2);
            
            if nStimTrials ~= nEEGTrials
            
                warning(['%s subject %d: %d stimulus trials but %d EEG trials. ' ...
                         'Skipping subject until trial mapping is verified.'], ...
                    datasetName,subjectIds(sub),nStimTrials,nEEGTrials);
            
                continue;
            
            end

            %% Ensuring that stim and eeg have the same length
            [stimFeature,eeg] = cndCutSameLength(stimFeature,eeg);

            %% TRF crossvalidation - determining optimal regularisation parameter
            % stimFeature.data = cellfun(@(x) x(1:eeg.fs*5,:),stimFeature.data,'UniformOutput',false);
            % eeg.data = cellfun(@(x) x(1:eeg.fs*5,:),eeg.data,'UniformOutput',false);
            if isfield(eeg,'origTrialPosition') && ~isempty(eeg.origTrialPosition) % if shuffled order
                trialIdx = find(eeg.origTrialPosition == ii);
            else
                trialIdx = ii;
            end
            if numel(trialIdx) ~= 1 || trialIdx > numel(eeg.data) || ...
                    trialIdx > numel(stimFeature.data)
                warning('%s subject %d: no unique matching trial %d; skipping.', ...
                    datasetName,subjectIds(sub),ii);
                continue
            end

            disp('Running mTRFcrossval')
            [stats,t] = mTRFcrossval(stimFeature.data(trialIdx), ...
                eeg.data(trialIdx),eeg.fs,dirTRF,tmin,tmax,lambdas, ...
                'split',3,'verbose',0);

            %% Calculating optimal lambda. Display and store results
            [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));
            disp(['r = ',num2str(maxR)])
            rAll(sub) = maxR;
            if isempty(rAllElec)
                rAllElec = nan(size(stats.r,3),nSubs);
            end
            rAllElec(:,sub) = squeeze(mean(stats.r(:,bestLambda,:),1));
            bestLambdaAll(sub) = lambdas(bestLambda);
            trialIdxAll(sub) = trialIdx;
            % origTrialPosition is retained so presentation order can be
            % checked before interpreting a trial-number adaptation curve.
            if isfield(eeg,'origTrialPosition') && ~isempty(eeg.origTrialPosition)
                presentationPositionAll(sub) = eeg.origTrialPosition(trialIdx);
            else
                presentationPositionAll(sub) = trialIdx; % assumed, not verified
            end

            %% Fit TRF model with optimal regularisation parameter
            % ORIGINAL LAB LOGIC: train using ALL trials, not only trialIdx.
            disp('Running mTRFtrain')
            model = mTRFtrain(stimFeature.data,eeg.data,eeg.fs,dirTRF, ...
                tmin,tmax,lambdas(bestLambda),'verbose',0);


            % TO REMOVE
            %% Check stimulus and EEG alignment before training
            
            nStimTrials = numel(stimFeature.data);
            nEEGTrials = numel(eeg.data);
            
            fprintf('Stimulus trials: %d | EEG trials: %d\n', ...
                nStimTrials,nEEGTrials);
            
            if nStimTrials ~= nEEGTrials
                error('Different numbers of stimulus and EEG trials.');
            end
            
            for tr = 1:nStimTrials
            
                nStimSamples = size(stimFeature.data{tr},1);
                nEEGSamples = size(eeg.data{tr},1);
            
                if nStimSamples ~= nEEGSamples
            
                    fprintf(['Mismatch in trial %d: ' ...
                        'STIM = %d samples, EEG = %d samples\n'], ...
                        tr,nStimSamples,nEEGSamples);
            
                end
            
            end
                
            % TO REMOVE


            %% Store TRF model
            nValid = nValid + 1;
            
            if nValid == 1
                modelAll = model;
            else
                modelAll(nValid) = model;
            end
            
            modelSubjectIds(nValid) = subjectIds(sub);
            plotEEG = eeg;
            downFs = eeg.fs;
            if isfield(eeg,'chanlocs') && ~isempty(eeg.chanlocs)
                chanlocs = eeg.chanlocs;
            end
            disp(['Mean r = ',num2str(mean(rAll(~isnan(rAll))))])
        end

        %% Save independent results for every dataset, feature and trial
        if nValid == 0
            warning('No valid subjects for %s trial %d.',datasetName,ii);
            continue
        end
        outputFilename = fullfile(outputDatasetFolder,sprintf('TRF_trial%03d.mat',ii));
        save(outputFilename,'chanlocs','modelAll','modelSubjectIds', ...
            'rAll','rAllElec','bestLambdaAll','trialIdxAll', ...
            'presentationPositionAll','subjectIds','bandpassFilterRange', ...
            'downFs','datasetName','nSubs','stimIdx','tmin','tmax', ...
            'lambdas','ii','preprocessedBand','-v7.3');
        fprintf('Saved: %s (%d valid subjects)\n',outputFilename,nValid);

        %% Optional original plotting (one figure per trial)
        if makePlots && ~isempty(plotEEG)
            figure('Position',[100,100,600,600]);
            % rAllElec columns are indexed by subject file. Select subjects
            % for which a model was actually trained.
            modelCols = find(~isnan(bestLambdaAll));
            if dirTRF == 1
                mTRF_plotForwardTRF(plotEEG,modelAll,rAllElec(:,modelCols));
            elseif dirTRF == -1
                mTRF_plotBackwardTRF(plotEEG,modelAll,rAllElec(:,modelCols));
            end
            saveas(gcf,fullfile(outputDatasetFolder,sprintf('TRF_trial%03d.png',ii)));
            close(gcf);
        end
    end
end

disp('Trial-by-trial envelope analysis completed.')
