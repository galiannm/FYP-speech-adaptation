% Cognition and Natural Sensory Processing (CNSP) Workshop
% Example 1 - Forward TRF
%
% This example script loads and preprocesses a publicly available dataset
% (you can use any of the dataset in the CNSP resources). Then, the script
% runs a typical forward TRF analysis.
%
% Note:
% This code was written with the assumption that all subjects were
% presented with the same set of stimuli. Hence, we use a single stimulus
% file (dataStim.mat) that applies to all subjects. This is compatible
% with scenarios with randomise presentation orders. In that case, the
% EEG/MEG trials should be sorted to match the single stimulus file. 
% The original order is preserved in a specific CND variable. If distinct
% subjects were presented with different stimuli, it is necessary to
% include a stimulus file per participant.
%
% CNSP-Workshop 2022
% https://cnspworkshop.net
% Author: Giovanni M. Di Liberto
% Copyright 2021 - Giovanni Di Liberto
%                  Nathaniel Zuk
%                  Michael Crosse
%                  Aaron Nidiffer
%                  Giorgia Cantisani
%                  (see license file for details)
% Last update: 24 June 2022

clear all
close all

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));

addpath(fullfile(projectRoot,'libs','cnsp_utils'));
addpath(fullfile(projectRoot,'libs','cnsp_utils','cnd'));
addpath(fullfile(projectRoot,'libs','mTRF-Toolbox_v2','mtrf'));
addpath(fullfile(projectRoot,'libs','NoiseTools'));
addpath(fullfile(projectRoot,'libs','eeglab'));

eeglab
% ======== ENVELOPE ========
% dataMainFolders = {'AAD_KULeuven/', 'AliceSpeech/', 'diliBach/', 'FDSpeech_L1/', 'FDSpeech_L2/', ...
%     'LalorNatSpeech/', 'LalorNatSpeechReverse/', 'SparrKULee1', 'childStoriesEEG_Sysoeva/', 'cocktailAttSwitch/', ...
%     'cocktailAttSwitch/', 'fNIRS_podcast/', 'musicImagery/', 'polyphonicBach/', 'spectrotemporalModulationsIgnacio/', ...
%     'hyperscanningAmir_Listening/', 'podcastListeningEEG/', 'ComedyShowEEG/', 'speechEmotionAdultEEG/', 'speechTrustEEG/', ....
%     'sungSpeechCantisani_melody', 'sungSpeechCantisani_song'};
% stimFeatureDataset = [1, 1, 1, 1, 1, ...
%                         1, 1, 1, 1, 1, ...
%                         2, 1, 1, 3, 1, ...
%                         1, 1, 1, 1, 1, ...
%                         1, 1];

dataMainFolders = {'AAD_KULeuven/', 'AliceSpeech/', 'diliBach/', 'FDSpeech_L1/', ...
    'LalorNatSpeech/', 'LalorNatSpeechReverse/', 'SparrKULee1', 'childStoriesEEG_Sysoeva/', 'cocktailAttSwitch/', ...
    'cocktailAttSwitch/', 'fNIRS_podcast/', 'musicImagery/', 'polyphonicBach/', 'spectrotemporalModulationsIgnacio/', ...
    'hyperscanningAmir_Listening/', 'podcastListeningEEG/', 'ComedyShowEEG/', 'speechEmotionAdultEEG/', 'speechTrustEEG/', ....
    'sungSpeechCantisani_melody', 'sungSpeechCantisani_song'};
stimFeatureDataset = [1, 1, 1, 1, ...
                        1, 1, 1, 1, 1, ...
                        2, 1, 1, 3, 1, ...
                        1, 1, 1, 1, 1, ...
                        1, 1];
feature = 'envelope';

% ======== WORD ONSET ========
% dataMainFolders = {'AliceSpeech/', 'FDSpeech_L1/', 'LalorNatSpeech/', 'childStoriesEEG_Sysoeva/', ...
%     'cocktailAttSwitch/', 'cocktailAttSwitch/', 'fNIRS_podcast/', 'hyperscanningAmir_Listening/', ...
%     'podcastListeningEEG/', 'ComedyShowEEG/', 'speechEmotionAdultEEG/', 'speechTrustEEG/'};
 
% stimFeatureDataset = [2, 5, 4, 2,  3, 4, 3, 4,  3, 11, 3, 3];

% feature = 'wordOnset';

toPreprocess = zeros(1,numel(dataMainFolders));


for iiDataset = 1:length(dataMainFolders)

    disp("*** "+dataMainFolders{iiDataset})
        
    %% Parameters - Natural speech listening experiment
    datasetName = char(erase(string(dataMainFolders{iiDataset}),"/"));
    dataMainFolder = fullfile(projectRoot,'datasets2',datasetName);
    dataCNDPath = fullfile(dataMainFolder,'dataCND');
    eegPreFolder = fullfile(dataCNDPath,'1-8Hz');

    reRefType = 'Avg'; % or 'Mastoids'
    bandpassFilterRange = [1,8]; % Hz (indicate 0 to avoid running the low-pass
                            % or high-pass filters or both)
                            % e.g., [0,8] will apply only a low-pass filter
                            % at 8 Hz
    downFs = 64; % Hz. *** fs/downFs must be an integer value ***

    eegFilenames = dir([eegPreFolder,'dataSub*.mat']);
    nSubs = length(eegFilenames);

    if downFs < bandpassFilterRange(2)*2
        disp('Warning: Be careful. The low-pass filter should use a cut-off frequency smaller than downFs/2')
    end

    %% Preprocess EEG - Natural speech listening experiment
    if toPreprocess(iiDataset)
    for sub = 1:nSubs
        % Loading EEG data
        eegFilename = [eegPreFolder,eegFilenames(sub).name];
        disp(['Loading EEG data: ',eegFilenames(sub).name])
        load(eegFilename,'eeg')
        
        eeg = cndNewOp(eeg,'Load'); % Saving the processing pipeline in the eeg struct
        
        % That dataset did not have a chanlocs in the CND. Fixing the issue
        if strcmp(dataMainFolder,'../datasets/spectrotemporalModulationsIgnacio/')
            load('chanlocs64.mat')
            eeg.chanlocs = chanlocs;
        end

        % Filtering - LPF (low-pass filter)
        if bandpassFilterRange(2) > 0
            hd = getLPFilt(eeg.fs,bandpassFilterRange(2));
            
            % A little coding trick - for loop vs cellfun
            if (0)
                % Filtering each trial/run with a for loop
                for ii = 1:length(eeg.data)
                    eeg.data{ii} = filtfilthd(hd,eeg.data{ii});
                end
            else
                % Filtering each trial/run with a cellfun statement
                eeg.data = cellfun(@(x) filtfilthd(hd,x),eeg.data,'UniformOutput',false);
            end
            
            % Filtering external channels
            if isfield(eeg,'extChan')
                for extIdx = 1:length(eeg.extChan)
                    eeg.extChan{extIdx}.data = cellfun(@(x) filtfilthd(hd,x),eeg.extChan{extIdx}.data,'UniformOutput',false);
                end
            end
            
            eeg = cndNewOp(eeg,'LPF');
        end
        
        % Downsampling EEG and external channels
        if downFs < eeg.fs
            eeg = cndDownsample(eeg,downFs);
        end
        
        % Filtering - HPF (high-pass filter)
        if bandpassFilterRange(1) > 0 
            hd = getHPFilt(eeg.fs,bandpassFilterRange(1));
            
            % Filtering EEG data
            eeg.data = cellfun(@(x) filtfilthd(hd,x),eeg.data,'UniformOutput',false);
            
            % Filtering external channels
            if isfield(eeg,'extChan')
                for extIdx = 1:length(eeg.extChan)
                    eeg.extChan{extIdx}.data = cellfun(@(x) filtfilthd(hd,x),eeg.extChan{extIdx}.data,'UniformOutput',false);
                end  
            end
            
            eeg = cndNewOp(eeg,'HPF');
        end
        
        % Replacing bad channels
        if isfield(eeg,'chanlocs')
            for tr = 1:length(eeg.data)
                eeg.data{tr} = removeBadChannels(eeg.data{tr}, eeg.chanlocs);
            end
        end
        
        % Re-referencing EEG data
        try
            eeg = cndReref(eeg,reRefType);
        catch
            reRefType = 'Avg';
            eeg = cndReref(eeg,'Avg');
        end
        
        % Removing initial padding (specific to this dataset)
        if isfield(eeg,'paddingStartSample')
            for tr = 1:length(eeg.data)
                eeg.data{tr} = eeg.data{tr}(1+eeg.paddingStartSample,:);
                for extIdx = 1:length(eeg.extChan)
                    eeg.extChan{extIdx}.data{tr} = eeg.extChan{extIdx}.data{tr}(1+eeg.paddingStartSample,:);
                end
            end
        end
        
        % Saving preprocessed data
        eegPreFilename = [dataMainFolder,dataCNDSubfolder,'pre_',eegFilenames(sub).name];
        disp(['Saving preprocessed EEG data: pre_',eegFilenames(sub).name])
        save(eegPreFilename,'eeg')
    end
    end

    %% The multivariate Temporal Response Function
    eegFilenames = dir(fullfile(eegPreFolder,'pre_dataSub*.mat'));
    eegFilenames = eegFilenames(~strcmp({eegFilenames.name},'pre_dataSub100.mat'));

    nSubs = numel(eegFilenames);

    if nSubs == 0
        warning('No preprocessed EEG found for %s; skipping.',datasetName);
        continue
    end

    subjectIds = nan(1,nSubs);
    for sub = 1:nSubs
        token = regexp(eegFilenames(sub).name, ...
            '^pre_dataSub(\d+)\.mat$','tokens','once');

        assert(~isempty(token), ...
            'Unexpected EEG filename: %s',eegFilenames(sub).name);

        subjectIds(sub) = str2double(token{1});
    end

    % TRF hyperparameters
    tmin = -100;
    tmax = 600;
    % lambdas = [1e-2,1e0,1e2]; % small set of lambdas (quick)
    lambdas = [1e-6,1e-3,1e-4,1e-3,1e-2,1e-1,1e0,1,1e2,1e3,1e4]; % larger set of lambdas (slower)
    dirTRF = 1; % Forward TRF model
    % Be careful: backward models (dirTRF-1) with many electrodes and large time-windows
    % can require long computational times. So, we suggest reducing the
    % dimensionality if you are just playing around with the code (e.g., select
    % only few electrodes and/or reduce the TRF window)
    stimIdx = stimFeatureDataset(iiDataset); 

    sharedStimFilename = fullfile(dataCNDPath,'dataStim.mat');
    subjectStimFiles = cell(1,nSubs);
    nStimTrials = zeros(1,nSubs);

    for sub = 1:nSubs
        if isfile(sharedStimFilename)
            subjectStimFiles{sub} = sharedStimFilename;
        else
            subjectStimFiles{sub} = fullfile(dataCNDPath, sprintf('dataStim%d.mat',subjectIds(sub)));
        end

        assert(isfile(subjectStimFiles{sub}), 'Missing stimulus file: %s',subjectStimFiles{sub});

        SubStim = load(subjectStimFiles{sub},'stim');

        assert(stimIdx <= size(SubStim.stim.data,1), 'Feature %d unavailable in %s',stimIdx,subjectStimFiles{sub});

        nStimTrials(sub) = size(SubStim.stim.data,2);

        if sub == 1
            fprintf('Dataset: %s | feature %d\n',datasetName,stimIdx);
            disp(SubStim.stim.names(stimIdx));
        end
    end

    nTrial = max(nStimTrials);
    if strcmp(feature,'envelope')
        resultDir = fullfile(scriptFolder,'resultsTrials_envelope');
    elseif strcmp(feature,'wordOnset')
        resultDir = fullfile(scriptFolder,'resultsTrials_wordOnset');
    else
        error('Unknown feature: %s',feature);
    end

    if ~isfolder(resultDir)
        mkdir(resultDir);
    end

    % TRF
    for ii = 1:nTrial
        modelAll = struct([]);
        modelSubjectIds = [];
        trialIdxAll = [];
        bestLambdaAll = [];
        rAll = [];
        rAllElec = [];
        chanlocs = [];
        subjectSamplingRates = [];

        % trialFig = figure('Name',sprintf('%s | trial %d',datasetName,ii), ...
        %     'Units','normalized','Position',[0.05,0.08,0.90,0.80]);
        
        for sub = 1:nSubs
            eegFilename = fullfile(eegFilenames(sub).folder,eegFilenames(sub).name);

            E = load(eegFilename,'eeg');
            eeg = E.eeg;

            Sstim = load(subjectStimFiles{sub},'stim');
            stim = Sstim.stim;

            % Use the EEG's actual sampling rate.
            eeg.fs = double(eeg.fs);
            stim.fs = double(stim.fs);

            % Selecting feature of interest ('stimIdx' feature)
            stimFeature = stim;
            % stimFeature.names = stimFeature.names{stimIdx};
            if iscell(stim.names)
                stimFeature.names = stim.names(stimIdx);
            elseif isstring(stim.names)
                stimFeature.names = stim.names(stimIdx);
            elseif ischar(stim.names)
                featureNames = cellstr(stim.names);
                assert(stimIdx <= numel(featureNames), ...
                    'Feature-name metadata does not match the selected feature.');
                stimFeature.names = featureNames(stimIdx);
            else
                error('Unsupported feature-name type: %s',class(stim.names));
            end
            stimFeature.data = stimFeature.data(stimIdx,:); % envelope or word onset

            if stimFeature.fs > eeg.fs
                stimFeature = cndDownsample(stimFeature,eeg.fs);
            end

            % Making sure that stim and neural data have the same length
            % The trial may end a few seconds after the end of the audio
            % e.g., the neural data may include the break between trials
            % It would be best to do this chunking at preprocessing, but let's
            % check here, just to be sure
            [stimFeature,eeg] = cndCheckStimNeural(stimFeature,eeg);
            
            % Standardise stim data (preserving the ratio between features)
            % This is thought for continuous signals e.g., speech envelope, eeg
            stimFeature = cndNormalise(stimFeature);
            % Standardise neural data (preserving the ratio between channels)
            eeg = cndNormalise(eeg);
                
            % Ensuring that stim and ieeg have the same length
            [stimFeature,eeg] = cndCutSameLength(stimFeature,eeg);

            % TRF crossvalidation - determining optimal regularisation parameter
            % stimFeature.data = cellfun(@(x) x(1:eeg.fs*5,:),stimFeature.data,'UniformOutput',false);
            % eeg.data = cellfun(@(x) x(1:eeg.fs*5,:),eeg.data,'UniformOutput',false);
            
            disp('Running mTRFcrossval')

            if isfield(eeg,'origTrialPosition') && ~isempty(eeg.origTrialPosition)
                trialIdx = find(eeg.origTrialPosition == ii);
            else
                if ii > numel(stimFeature.data) || ii > numel(eeg.data)
                    warning('%s: subject %d has no trial %d; skipping.', ...
                        datasetName,subjectIds(sub),ii);
                    continue
                end

                trialIdx = ii;
            end

            if isempty(trialIdx)
                warning('%s: subject %d has no presentation position %d; skipping.', ...
                    datasetName,subjectIds(sub),ii);
                continue
            end

            assert(numel(trialIdx) == 1, ...
                '%s: subject %d, position %d has %d matching trials.', ...
                datasetName,subjectIds(sub),ii,numel(trialIdx));

            assert(trialIdx >= 1 && trialIdx <= numel(stimFeature.data), ...
                'Selected trial is outside the available stimulus trials.');

            [stats,t] = mTRFcrossval(stimFeature.data(trialIdx),eeg.data(trialIdx),eeg.fs,dirTRF,tmin,tmax,lambdas,'split',3,'verbose',0);

            % Calculating optimal lambda. Display and store results
            [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));
            disp(['r = ',num2str(maxR)]);
            
            % Fit TRF model with optimal regularisation parameter
            disp('Running mTRFtrain');
            model = mTRFtrain(stimFeature.data,eeg.data,eeg.fs,dirTRF,tmin,tmax,lambdas(bestLambda),'verbose',0);

            if isempty(modelAll)
                modelAll = model;
            else
                modelAll(end+1) = model;
            end
            modelSubjectIds(end+1) = subjectIds(sub);
            trialIdxAll(end+1) = trialIdx;
            bestLambdaAll(end+1) = lambdas(bestLambda);
            subjectSamplingRates(end+1) = eeg.fs;

            rAll(end+1) = maxR;
            rAllElec(:,end+1) = squeeze(mean(stats.r(:,bestLambda,:),1));

            if isfield(eeg,'chanlocs') && ~isempty(eeg.chanlocs)
                chanlocs = eeg.chanlocs;
            else
                chanlocs = [];
            end

            % figure(trialFig);
            % if dirTRF == 1
            %     mTRF_plotForwardTRF(eeg,modelAll,rAllElec);
            % elseif dirTRF == -1
            %     mTRF_plotBackwardTRF(eeg,modelAll,rAllElec);
            % end
            % drawnow;
        end
        if isempty(modelAll)
            warning('No valid subjects for %s trial %d.',datasetName,ii);
            close(trialFig);
            continue
        end
        outputFilename = fullfile(resultDir,sprintf( ...
            'trfOut_Trial_%s_stim%d_trial%03d.mat', ...
            datasetName,stimIdx,ii));

        save(outputFilename, ...
            'chanlocs','modelAll','modelSubjectIds', ...
            'rAll','rAllElec','trialIdxAll','bestLambdaAll', ...
            'subjectSamplingRates','bandpassFilterRange', ...
            'reRefType','dataMainFolder','datasetName', ...
            'nSubs','stimIdx','tmin','tmax','lambdas','ii');

        fprintf('Saved: %s\n',outputFilename);
    end
end 

