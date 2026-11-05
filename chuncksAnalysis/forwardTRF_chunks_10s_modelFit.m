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

addpath ../../libs/cnsp_utils
addpath ../../libs/cnsp_utils/cnd
addpath ../../libs/mTRF-Toolbox_v2/mtrf
addpath ../../libs/NoiseTools
addpath ../../libs/eeglab
eeglab

% envelope
% dataMainFolders = {'AAD_KULeuven', 'AliceSpeech', 'diliBach/', 'FDSpeech_L1', 'FDSpeech_L2', ...
%     'LalorNatSpeech/', 'LalorNatSpeechReverse/', 'SparrKULee1', 'childStoriesEEG_Sysoeva', 'cocktailAttSwitch', ...
%     'cocktailAttSwitch/', 'fNIRS_podcast', 'musicImagery', 'polyphonicBach/', 'spectrotemporalModulationsIgnacio/'};

% word onset
dataMainFolders = {'AliceSpeech', 'FDSpeech_L1', 'LalorNatSpeech', 'childStoriesEEG_Sysoeva', ...
'cocktailAttSwitch', 'cocktailAttSwitch', 'fNIRS_podcast', 'hyperscanningAmir_Listening', ...
'podcastListeningEEG', 'ComedyShowEEG', 'speechEmotionAdultEEG', 'speechTrustEEG', };

% dataMainFolders = {'LalorNatSpeech/','LalorNatSpeechReverse/'};
% dataMainFolders = {'LalorNatSpeech/','LalorNatSpeechReverse/','diliBach/','emilyDialogue/','AliceSpeech/','cocktailAttSwitch/','anxietyEEG/'};
% dataMainFolders = {'LalorNatSpeech/','LalorNatSpeechReverse/','diliBach/','laughterEEG/','cocktailAttSwitch/','cocktailAttSwitch/','spectrotemporalModulationsIgnacio/','polyphonicBach/','polyphonicBach/'};
% dataMainFoldersLabels = {'LalorNatSpeech','LalorNatSpeechReverse','diliBach','laughterEEG','cocktailAttSwitch_target','cocktailAttSwitch_masker','spectrotempMod-clean','spectrotempMod-temp60','spectrotempMod-temp80','spectrotempMod-spectr60','spectrotempMod-spectr80','polyphonicBachMain/','polyphonicBachSupport/'};
% dataMainFolders = {'spectrotemporalModulationsIgnacio/'};
% toPreprocess = [0];
% dataMainFolders = {'musicImagery/'};
toPreprocess = zeros(1, length(dataMainFolders));

% envelope
% stimFeatureDataset = [1,1,1,1,1, ...
%                         1,1,1,1,1, ...
%                         2,1,1,3,1];

% word onset
stimFeatureDataset = [2, 5, 4, 2, 3, 4, 3, 4, 3, 11, 3, 3];

% stimFeatureDataset = [2,2];
% dataMainFolders = {'laughterEEG/','musicImagery/'};
% toPreprocess = [0,0];
for iiDataset = 1:length(dataMainFolders)

disp("*** "+dataMainFolders{iiDataset})
    
%% Parameters - Natural speech listening experiment
scriptFolder = fileparts(mfilename('fullpath'));
dataMainFolder = fullfile( ...
    scriptFolder, '..', 'datasets2', dataMainFolders{iiDataset});
dataCNDSubfolder = 'dataCND/';
resultsDir = fullfile(scriptFolder, 'resultsChunks', 'wordOnset');

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

reRefType = 'Avg'; % or 'Mastoids'
bandpassFilterRange = [1,8]; % Hz (indicate 0 to avoid running the low-pass
                          % or high-pass filters or both)
                          % e.g., [0,8] will apply only a low-pass filter
                          % at 8 Hz
% if iiDataset >= 8
%     downFs = 125; % Hz. *** fs/downFs must be an integer value ***
% else
%     downFs = 64;
% end

downFs = [];

eegFilenames = dir(fullfile(dataMainFolder, ...
    dataCNDSubfolder, '1-8Hz', 'pre_dataSub*.mat'));
nSubs = length(eegFilenames);

if downFs < bandpassFilterRange(2)*2
    disp('Warning: Be careful. The low-pass filter should use a cut-off frequency smaller than downFs/2')
end

%% Preprocess EEG - Natural speech listening experiment
if toPreprocess(iiDataset)
   
for sub = 1:nSubs
    % Loading EEG data
    eegFilename = [dataMainFolder,dataCNDSubfolder,eegFilenames(sub).name];
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
    % if downFs < eeg.fs
    %     eeg = cndDownsample(eeg,downFs);
    % end
    
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
    eegPreFilename = fullfile(dataMainFolder, dataCNDSubfolder, '1-8Hz', eegFilenames(sub).name);   
    disp(['Saving preprocessed EEG data: pre_',eegFilenames(sub).name])
    save(eegPreFilename,'eeg')
end
end

%% The multivariate Temporal Response Function
% 0-10s, 10-20s, etc
% Loading Stimulus data
% stimFilename = [dataMainFolder,dataCNDSubfolder,'dataStim.mat'];
% if exist(stimFilename, 'file')==1
%     disp(['Loading stimulus data: ','dataStim.mat'])
%     singleDataStim = true;
% else
%     stimFilename = [dataMainFolder,dataCNDSubfolder,'dataStim1.mat'];
%     disp(['Loading stimulus data: ','dataStim1.mat'])
%     singleDataStim = false;
% end
stimFolder = fullfile(dataMainFolder, dataCNDSubfolder);
sharedStimFile = fullfile(stimFolder, 'dataStim.mat');
if isfile(sharedStimFile)
    singleDataStim = true;
    stimFilename = sharedStimFile;
else
    singleDataStim = false;
    stimFiles = dir(fullfile(stimFolder, 'dataStim*.mat'));

    if isempty(stimFiles)
        error('No stimulus files found in %s', stimFolder);
    end
    
    stimFilename = fullfile(stimFolder, stimFiles(1).name);
end
load(stimFilename,'stim')

% Downsampling stim if necessary
stim.data = cellfun(@double, stim.data, 'UniformOutput', false);
stim.fs = double(stim.fs);
downFs = double(downFs);
% if downFs < stim.fs
%     stim = cndDownsample(stim,downFs);
% end

if ~isfield(stim,'condIdxs')
    stim.condIdxs = ones(size(stim.data,2),1);
end

for cond2select = 1:length(unique(stim.condIdxs))
for ii = 1:20
    ii
% TRF hyperparameters
tmin = -100;
% tmax = 400;
tmax = 700;
% if iiDataset == 7
%     tmin = tmin-500;
%     tmax = tmax-500;
% end
% lambdas = [1e-2,1e0,1e2]; % small set of lambdas (quick)
lambdas = [1e-6,1e-3,1e-4,1e-3,1e-2,1e-1,1e0,1,1e2,1e3,1e4]; % larger set of lambdas (slower)
dirTRF = 1; % Forward TRF model
% Be careful: backward models (dirTRF-1) with many electrodes and large time-windows
% can require long computational times. So, we suggest reducing the
% dimensionality if you are just playing around with the code (e.g., select
% only few electrodes and/or reduce the TRF window)
stimIdx = stimFeatureDataset(iiDataset); 

if singleDataStim
    % Loading Stimulus data
    % stimFilename = [dataMainFolder,dataCNDSubfolder,'dataStim.mat'];
    % disp(['Loading stimulus data: ','dataStim.mat'])
    % load(stimFilename,'stim')
    load(stimFilename, 'stim');
    stim.fs = double(stim.fs);
    % Downsampling stim if necessary
    % if downFs < stim.fs
    %     stim = cndDownsample(stim,downFs);
    % end
    if ~isfield(stim,'condIdxs')
        stim.condIdxs = ones(size(stim.data,2),1);
    end
end

% TRF
rAll = [];
rAllElec = [];
modelAll = struct([]);
includedSubjectIds = [];
% figure('Position',[100,100,600,600]);
for sub = 1:nSubs
    % Loading preprocessed EEG
    subIdx = regexp(eegFilenames(sub).name, '\d+', 'match');
    subIdx = str2double(subIdx{1});
    eegPreFilename = fullfile(dataMainFolder, dataCNDSubfolder, '1-8Hz', eegFilenames(sub).name);    
    disp(['Loading preprocessed EEG data: ',eegFilenames(sub).name])
    load(eegPreFilename,'eeg')

    if ~singleDataStim
        % Loading Stimulus data
        stimFilename = fullfile(dataMainFolder, dataCNDSubfolder, ['dataStim', num2str(subIdx), '.mat']);
        disp(['Loading stimulus data: ','dataStim',num2str(subIdx),'.mat'])
        load(stimFilename,'stim')
        stim.fs = double(stim.fs);
        % Downsampling stim if necessary
        % if downFs < stim.fs
        %     stim = cndDownsample(stim,downFs);
        % end
    end
    

    % Selecting feature of interest ('stimIdx' feature)
    % stimFeature = stim;
    % stimFeature.names = stimFeature.names{stimIdx};
    % stimFeature.data = stimFeature.data(stimIdx,:); % envelope or word onset

    % Selecting feature of interest ('stimIdx' feature)
    stimFeature = stim;
    % stimFeature.names = stimFeature.names{stimIdx};
    if iscell(stimFeature.names)
        stimFeature.names = stimFeature.names{stimIdx};
    else
        stimFeature.names = stimFeature.names;
    end
    stimFeature.data = stimFeature.data(stimIdx,:); % envelope or word onset
    if ~isfield(stimFeature, 'condIdxs') || isempty(stimFeature.condIdxs)
        stimFeature.condIdxs = ones(1, numel(stimFeature.data));
    end
    
    % Match stimulus sampling frequency to preprocessed EEG
    eeg.fs = double(eeg.fs);
    stimFeature.fs = double(stimFeature.fs);
    
    if stimFeature.fs ~= eeg.fs
    
        fprintf('Resampling STIM from %g Hz to %g Hz\n', ...
            stimFeature.fs, eeg.fs);
    
        oldFs = stimFeature.fs;
        newFs = eeg.fs;
    
        for tr = 1:numel(stimFeature.data)
            stimFeature.data{tr} = resample( ...
                double(stimFeature.data{tr}), newFs, oldFs);
        end
    
        stimFeature.fs = newFs;
    
    end
    
    assert(stimFeature.fs == eeg.fs, ...
        'EEG and STIM sampling frequencies do not match');
    
    % Save actual sampling frequency
    downFs = eeg.fs;
    
    fprintf('EEG fs = %g Hz | STIM fs = %g Hz\n', ...
        eeg.fs, stimFeature.fs);
    
    % Making sure that stim and neural data have the same length
    % The trial may end a few seconds after the end of the audio
    % e.g., the neural data may include the break between trials
    % It would be best to do this chunking at preprocessing, but let's
    % check here, just to be sure
    % Skip subjects where EEG and stimulus trial counts do not match
    if numel(stimFeature.data) ~= numel(eeg.data)
        warning('Skipping subject %d: STIM has %d trials, EEG has %d trials', ...
            subIdx, numel(stimFeature.data), numel(eeg.data));
        continue
    end
    [stimFeature,eeg] = cndCheckStimNeural(stimFeature,eeg);
    
    % Standardise stim data (preserving the ratio between features)
    % This is thought for continuous signals e.g., speech envelope, eeg
    stimFeature = cndNormalise(stimFeature);
    % Standardise neural data (preserving the ratio between channels)
    eeg = cndNormalise(eeg);
        
    % Ensuring that stim and ieeg have the same length
    [stimFeature,eeg] = cndCutSameLength(stimFeature,eeg);
    
    % TRF crossvalidation - determining optimal regularisation parameter
%     stimFeature.data = cellfun(@(x) x(1:eeg.fs*5,:),stimFeature.data,'UniformOutput',false);
%     eeg.data = cellfun(@(x) x(1:eeg.fs*5,:),eeg.data,'UniformOutput',false);
    sampStart = 1+eeg.fs*10*(ii-1);
    sampEnd = eeg.fs*10*(ii);
    toBreak = false; 
    for jjTr = 1:length(eeg.data)
        if sampEnd > size(eeg.data{jjTr},1)
            toBreak = true;
        end
    end
    if toBreak
        disp('Reached the end of the shortest trial')
        break
    end
    
    stimFeature.data = cellfun(@(x) x(sampStart:sampEnd,:),stimFeature.data,'UniformOutput',false);
    eeg.data = cellfun(@(x) x(sampStart:sampEnd,:),eeg.data,'UniformOutput',false);

    disp('Running mTRFcrossval')
    conditions = unique(stimFeature.condIdxs, 'stable');

    % Skip subjects without the requested condition
    if cond2select > numel(conditions)
        warning('Skipping subject %d: condition %d unavailable', ...
            subIdx, cond2select);
        continue
    end
    
    if iscell(stimFeature.condIdxs)
        trials2keep = strcmp(stimFeature.condIdxs, conditions{cond2select});
    else
        trials2keep = stimFeature.condIdxs == conditions(cond2select);
    end
    % [stats,t] = mTRFcrossval(stimFeature.data(stimFeature.condIdxs==cond2select),eeg.data(stimFeature.condIdxs==cond2select),eeg.fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
    [stats,t] = mTRFcrossval( ...
        stimFeature.data(trials2keep), ...
        eeg.data(trials2keep), ...
        eeg.fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
    % Calculating optimal lambda. Display and store results
    [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));
    disp(['r = ',num2str(maxR)])
    % rAll(sub) = maxR;
    % rAllElec(:,sub) = squeeze(mean(stats.r(:,bestLambda,:),1));
    rAll(end+1) = maxR;
    rAllElec(:,end+1) = squeeze(mean(stats.r(:,bestLambda,:),1));
    
    % Fit TRF model with optimal regularisation parameter
    disp('Running mTRFtrain')
    % model = mTRFtrain(stimFeature.data,eeg.data,eeg.fs,dirTRF,tmin,tmax,lambdas(bestLambda),'verbose',0);
    
    conditions = unique(stimFeature.condIdxs, 'stable');

    if iscell(stimFeature.condIdxs)
        trials2keep = strcmp(stimFeature.condIdxs, conditions{cond2select});
    else
        trials2keep = stimFeature.condIdxs == conditions(cond2select);
    end
    model = mTRFtrain( ...
        stimFeature.data(trials2keep), ...
        eeg.data(trials2keep), ...
        eeg.fs,dirTRF,tmin,tmax,lambdas(bestLambda),'verbose',0);
    
    % Store TRF model
    if isempty(modelAll)
        modelAll = model;
    else
        modelAll(end+1) = model;
    end
    includedSubjectIds(end+1) = subIdx;
    
    % if dirTRF == 1
    %     mTRF_plotForwardTRF(eeg,modelAll,rAllElec);
    % elseif dirTRF == -1
    %     mTRF_plotBackwardTRF(eeg,modelAll,rAllElec);
    % end
    
    disp(['Mean r = ',num2str(mean(rAll))])
    
    % drawnow;
end

if isfield(eeg,'chanlocs') && ~isempty(eeg.chanlocs)
    chanlocs = eeg.chanlocs;
else
    chanlocs = [];
end
if ~isempty(modelAll)

    filename = sprintf( ...
        'trfOut_10sChunk_%s_cond%d_stim%d_%d.mat', ...
        char(strip(string(dataMainFolders{iiDataset}), '/')), ...
        cond2select, stimIdx, ii);

    save(fullfile(resultsDir, filename), ...
        'chanlocs', 'modelAll', 'rAllElec', ...
        'bandpassFilterRange', 'downFs', 'reRefType', ...
        'dataMainFolder', 'nSubs', 'stimIdx', ...
        'tmin', 'tmax', 'lambdas');

end
if toBreak
    break
end
end
end
end
