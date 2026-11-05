%% TEST FORWARD TRF FITTING - standalone, does not change the main script
% Save this file in CNSP/myScripts/AnalysisScripts/
% Tests AliceSpeech, LalorNatSpeech and cocktailAttSwitch one by one.
% Uses only the FIRST 10-second chunk and keeps each dataset's results separate.

clear;
close all;

%% Paths (independent of MATLAB's current working directory)
scriptFolder = fileparts(mfilename('fullpath'));
rootFolder = fullfile(scriptFolder, '..', '..');
libsFolder = fullfile(rootFolder, 'libs');

addpath(fullfile(libsFolder, 'cnsp_utils'));
addpath(fullfile(libsFolder, 'cnsp_utils', 'cnd'));
addpath(fullfile(libsFolder, 'mTRF-Toolbox_v2', 'mtrf'));
addpath(fullfile(libsFolder, 'NoiseTools'));
addpath(fullfile(libsFolder, 'eeglab'));
% EEGLAB graphical interface is not required for fitting.
% Uncomment if a CNSP helper needs EEGLAB initialization:
% eeglab;

%% Test configuration (no changes needed between datasets)
tests = {
    'speechEmotionAdultEEG',       3;  % Word Onset 
    'speechTrustEEG',    3;  % Word Onset 
};

featureName = 'wordOnset';
chunkLength = 10;  % seconds
maxChunks = 1;    % quick smoke test: only 0-10 s
maxSubjects = Inf; % all available participants; set to 2 for an even faster check

tmin = -100;
tmax = 600;
lambdas = [1e-6,1e-3,1e-4,1e-2,1e-1,1e0,1e2,1e3,1e4];
dirTRF = 1;
bandpassFilterRange = [1 8];

% Separate from the main resultsChunks folder: NEVER overwrites main results.
resultDir = fullfile(scriptFolder, 'testResultsChunks', featureName);
if ~exist(resultDir, 'dir'), mkdir(resultDir); end

summaryDataset = strings(size(tests,1),1);
summaryStatus = strings(size(tests,1),1);
summarySubjects = zeros(size(tests,1),1);
summaryDetail = strings(size(tests,1),1);

%% Run tests separately; a failure in one dataset does not stop the others
for iTest = 1:size(tests,1)
    datasetName = tests{iTest,1};
    stimIdx = tests{iTest,2};
    summaryDataset(iTest) = string(datasetName);
    summaryStatus(iTest) = "FAILED";

    fprintf('\n========== TEST %d/%d: %s (stimIdx=%d) ==========\n', ...
        iTest, size(tests,1), datasetName, stimIdx);

    try
        dataCNDfolder = fullfile(rootFolder, 'datasets2', datasetName, 'dataCND');
        eegFolder = fullfile(dataCNDfolder, '1-8Hz');
        eegFilenames = dir(fullfile(eegFolder, 'pre_dataSub*.mat'));
        if isempty(eegFilenames)
            error('No pre_dataSub*.mat in %s', eegFolder);
        end
        fprintf('Found %d subjects\n', numel(eegFilenames));

        sharedStimFile = fullfile(dataCNDfolder, 'dataStim.mat');
        singleDataStim = isfile(sharedStimFile);
        if singleDataStim
            loaded = load(sharedStimFile, 'stim');
            sharedStim = loaded.stim;
            fprintf('Using shared stimulus dataStim.mat\n');
        else
            sharedStim = []; % subject-specific dataStimN.mat files
            fprintf('Using subject-specific stimulus files\n');
        end

        for ii = 1:maxChunks
            rAll = [];
            rAllElec = [];
            modelAll = struct([]);
            includedSubjectIdx = [];
            wordCountMeanAll = [];
            wordCountStdAll = [];
            chanlocs = [];
            downFs = NaN;
            reRefType = 'not applied here; inspect EEG preprocessingPipeline';
            selectedCondition = 1; % single-condition tests only

            nToTest = min(numel(eegFilenames), maxSubjects);
            for sub = 1:nToTest
                try
                    eegName = eegFilenames(sub).name;
                    idMatch = regexp(eegName, '^pre_dataSub(\d+)\.mat$', 'tokens', 'once');
                    if isempty(idMatch), error('Unexpected filename: %s', eegName); end
                    subIdx = str2double(idMatch{1});
                    fprintf('\n%s | chunk %d | subject %d\n', datasetName, ii, subIdx);

                    loaded = load(fullfile(eegFolder, eegName), 'eeg');
                    eeg = loaded.eeg;
                    eeg.fs = double(eeg.fs);

                    if singleDataStim
                        stim = sharedStim;
                    else
                        stimFilename = fullfile(dataCNDfolder, sprintf('dataStim%d.mat', subIdx));
                        if ~isfile(stimFilename)
                            warning('Missing stimulus file for subject %d; skipping', subIdx);
                            continue;
                        end
                        loaded = load(stimFilename, 'stim');
                        stim = loaded.stim;
                    end
                    stim.fs = double(stim.fs);
                    if stimIdx > size(stim.data, 1)
                        error('Feature index %d not found (%d features)', stimIdx, size(stim.data, 1));
                    end

                    % Extract exactly one feature, just as in the main script.
                    stimFeature = stim;
                    if isfield(stimFeature,'names') && iscell(stimFeature.names)
                        stimFeature.names = stimFeature.names{stimIdx};
                    end
                    stimFeature.data = stimFeature.data(stimIdx,:);

                    % Condition labels: missing means all trials are condition 1.
                    if ~isfield(stimFeature,'condIdxs') || isempty(stimFeature.condIdxs)
                        stimFeature.condIdxs = ones(1,numel(stimFeature.data));
                    elseif iscell(stimFeature.condIdxs)
                        if all(cellfun(@(x) isnumeric(x) && isscalar(x),stimFeature.condIdxs))
                            stimFeature.condIdxs = cell2mat(stimFeature.condIdxs);
                        else
                            error('Non-numeric cell condition labels: inspect manually');
                        end
                    end
                    stimFeature.condIdxs = stimFeature.condIdxs(:)';
                    if numel(stimFeature.condIdxs) ~= numel(stimFeature.data)
                        error('Mismatch between stimulus trial count and condition labels');
                    end

                    % For these three smoke tests the desired condition is 1.
                    trials2keep = stimFeature.condIdxs == selectedCondition;
                    if nnz(trials2keep) < 2
                        warning('Only %d trials in condition %d; skipping', ...
                            nnz(trials2keep), selectedCondition);
                        continue;
                    end

                    % Word counts BEFORE resampling or normalization.
                    originalFs = double(stimFeature.fs);
                    wordStart = 1 + round(originalFs*chunkLength*(ii-1));
                    wordEnd = round(originalFs*chunkLength*ii);
                    wordCounts = nan(1,numel(stimFeature.data));
                    for tr = find(trials2keep)
                        wordSignal = double(stimFeature.data{tr});
                        if wordEnd <= size(wordSignal,1)
                            wordCounts(tr) = nnz(wordSignal(wordStart:wordEnd,:));
                        end
                    end
                    fprintf('Selected feature size (trial 1): %s\n', ...
                        mat2str(size(stimFeature.data{1})));

                    if numel(stimFeature.data) ~= numel(eeg.data)
                        warning('STIM has %d trials, EEG has %d; skipping', ...
                            numel(stimFeature.data),numel(eeg.data));
                        continue;
                    end

                    % Resample stimulus to existing preprocessed EEG rate.
                    if stimFeature.fs ~= eeg.fs
                        fprintf('Resampling STIM %g -> EEG %g Hz\n',stimFeature.fs,eeg.fs);
                        for tr=1:numel(stimFeature.data)
                            stimFeature.data{tr} = resample(double(stimFeature.data{tr}), ...
                                eeg.fs, double(stimFeature.fs));
                        end
                        stimFeature.fs = eeg.fs;
                    end

                    [stimFeature,eeg] = cndCheckStimNeural(stimFeature,eeg);
                    stimFeature = cndNormalise(stimFeature);
                    eeg = cndNormalise(eeg);
                    [stimFeature,eeg] = cndCutSameLength(stimFeature,eeg);

                    sampStart = 1 + round(eeg.fs*chunkLength*(ii-1));
                    sampEnd = round(eeg.fs*chunkLength*ii);
                    if any(cellfun(@(x) size(x,1) < sampEnd, eeg.data(trials2keep))) || ...
                       any(cellfun(@(x) size(x,1) < sampEnd, stimFeature.data(trials2keep)))
                        warning('Selected condition contains short trials; skipping subject');
                        continue;
                    end
                    stimFeature.data = cellfun(@(x) x(sampStart:sampEnd,:), ...
                        stimFeature.data(trials2keep),'UniformOutput',false);
                    eeg.data = cellfun(@(x) x(sampStart:sampEnd,:), ...
                        eeg.data(trials2keep),'UniformOutput',false);

                    fprintf('Running mTRFcrossval...\n');
                    [stats,~] = mTRFcrossval(stimFeature.data,eeg.data, ...
                        eeg.fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
                    lambdaScores = squeeze(mean(mean(stats.r,1,'omitnan'),3,'omitnan'));
                    [maxR,bestLambda] = max(lambdaScores);
                    fprintf('r = %.5f | best lambda = %g\n',maxR,lambdas(bestLambda));

                    fprintf('Running mTRFtrain...\n');
                    model = mTRFtrain(stimFeature.data,eeg.data, ...
                        eeg.fs,dirTRF,tmin,tmax,lambdas(bestLambda),'verbose',0);

                    if isempty(modelAll)
                        modelAll = model;
                        if isfield(eeg,'chanlocs'), chanlocs=eeg.chanlocs; end
                    else
                        if ~isequal(size(model.w),size(modelAll(1).w))
                            warning('Subject %d has incompatible model dimensions; skipping',subIdx);
                            continue;
                        end
                        modelAll(end+1) = model;
                    end
                    rAll(end+1) = maxR;
                    rAllElec(:,end+1) = squeeze(mean(stats.r(:,bestLambda,:),1,'omitnan'));
                    includedSubjectIdx(end+1) = subIdx;
                    selectedWordCounts = wordCounts(trials2keep);
                    wordCountMeanAll(end+1) = mean(selectedWordCounts,'omitnan');
                    wordCountStdAll(end+1) = std(selectedWordCounts,0,'omitnan');
                    downFs = eeg.fs;
                    fprintf('Words/trial: %.2f +/- %.2f\n', ...
                        wordCountMeanAll(end),wordCountStdAll(end));
                catch subjectErr
                    warning('%s, subject %d: %s',datasetName,sub,subjectErr.message);
                end
            end

            if isempty(modelAll)
                fprintf('\n%s: NO MODELS FITTED for chunk %d\n',datasetName,ii);
                continue;
            end

            % Test output ONLY: does not touch resultsChunks/wordOnset.
            filename = sprintf('TEST_trfOut_%dsChunk_%s_cond%d_stim%d_%d.mat', ...
                chunkLength,datasetName,selectedCondition,stimIdx,ii);
            outputFile = fullfile(resultDir,filename);
            nSubs = numel(eegFilenames);
            dataMainFolder = fullfile(rootFolder,'datasets2',datasetName);
            save(outputFile,'chanlocs','modelAll','rAll','rAllElec', ...
                'bandpassFilterRange','downFs','reRefType', ...
                'dataMainFolder','nSubs','stimIdx','tmin','tmax','lambdas', ...
                'includedSubjectIdx','wordCountMeanAll','wordCountStdAll', ...
                'selectedCondition','chunkLength');

            % Reload from disk: make sure file really contains expected results.
            s = load(outputFile);
            assert(numel(s.modelAll)==numel(s.includedSubjectIdx));
            assert(numel(s.modelAll)==numel(s.wordCountMeanAll));
            assert(numel(s.modelAll)==numel(s.wordCountStdAll));
            fprintf('\nPASS: %s | fitted %d/%d participants | %s\n', ...
                datasetName,numel(s.modelAll),nToTest,outputFile);
            summaryStatus(iTest) = "PASS";
            summarySubjects(iTest) = numel(s.modelAll);
            summaryDetail(iTest) = string(outputFile);
        end
    catch datasetErr
        summaryDetail(iTest) = string(datasetErr.message);
        warning('Dataset test failed for %s: %s',datasetName,datasetErr.message);
    end
end

%% Summary
fprintf('\n========== TEST SUMMARY ==========\n');
disp(table(summaryDataset,summaryStatus,summarySubjects,summaryDetail, ...
    'VariableNames',{'Dataset','Status','ModelsFitted','Details'}));
