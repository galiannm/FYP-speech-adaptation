
%% Analyse results

clear all
close all

addpath ../libs/cnsp_utils
addpath ../libs/cnsp_utils/cnd
addpath ../libs/mTRF-Toolbox_v2/mtrf
addpath ../libs/NoiseTools
addpath ../libs/eeglab
eeglab

%%
defaultChunkLength = 100; % Seconds
close all
cmaxPred = 0.1;

superSub = true;

% Path to saved results
scriptFolder = fileparts(mfilename('fullpath'));

if superSub 
    resultsDir = fullfile(scriptFolder, 'resultsChunks', 'superSubject', 'wordOnset');
else
    resultsDir = fullfile(scriptFolder, 'resultsChunks', 'wordOnset');
end

% Dataset names matching the saved results
% Discover datasets automatically because word-onset feature indices
% can differ across datasets.
resultFiles = dir(fullfile(resultsDir, 'trfOut_100sChunk_*.mat')); 

dataMainFolders = {};

for iFile = 1:length(resultFiles)

    tokens = regexp(resultFiles(iFile).name, ...
        '^trfOut_100sChunk_(.*)_(\d+)\.mat$', 'tokens', 'once');

    if isempty(tokens)
        continue
    end

    datasetName = [tokens{1}, '_'];

    if ~ismember(datasetName, dataMainFolders)
        dataMainFolders{end+1} = datasetName;
    end
end

dataMainFoldersLabels = dataMainFolders;

if isempty(dataMainFolders)
    error('No word-onset result files found in %s', resultsDir);
end

nChunks = 20;
clear n400Val n400ValRaw
clear dataLen

% Initialise results
nDatasets = length(dataMainFolders);

n400Val = nan(nDatasets,nChunks);
n400ValRaw = nan(nDatasets,nChunks);

n400Peak = nan(nDatasets,nChunks);
n400PeakRaw = nan(nDatasets,nChunks);

predcorr = nan(nDatasets,nChunks);
pcorr_nonnorm = nan(nDatasets,nChunks);

pcorr_posterior = nan(nDatasets,nChunks);

dataLen = zeros(nDatasets,1);

% N400 latency window
n400Start = 300;
n400End = 500;

for iiDataset = 1:length(dataMainFolders)

    fprintf('\nProcessing %s\n', ...
        dataMainFoldersLabels{iiDataset});

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks

        filename = "trfOut_100sChunk_" + ...
            dataMainFolders{iiDataset} + ii + ".mat";

        filepath = fullfile(resultsDir,filename);

        if ~isfile(filepath)
            fprintf('Missing chunk: %s\n',filename);
            break
        end

        S = load(filepath, ...
            'chanlocs','modelAll','rAllElec', ...
            'bandpassFilterRange','downFs','reRefType', ...
            'dataMainFolder','nSubs','stimIdx', ...
            'tmin','tmax','lambdas');

        if ~isfield(S,'modelAll') || ~isfield(S,'rAllElec')
            warning('Missing model or correlations: %s',filename);
            break
        end

        modelAll = S.modelAll;
        rAllElec = S.rAllElec;

        if isempty(modelAll) || isempty(rAllElec)
            fprintf('Empty models: %s\n',filename);
            break
        end

        % Plot average TRF
        plotNormFlag = 0;

        isThereData = zeros(length(modelAll),1);

        for iiModel = 1:length(modelAll)
            isThereData(iiModel) = ...
                ~isempty(modelAll(iiModel).w);
        end

        modelAll = modelAll(logical(isThereData));

        if isempty(modelAll)
            fprintf('No valid models: %s\n',filename);
            break
        end

        avgModel = mTRFmodelAvg(modelAll,plotNormFlag);
        avgModel.w = squeeze(avgModel.w);

        t = avgModel.t(:);
        w = avgModel.w;

        % Ensure TRF weights are latency x electrode
        if size(w,1) ~= numel(t)
            if size(w,2) == numel(t)
                w = w';
            else
                warning('Unexpected TRF dimensions: %s', filename);
                break
            end
        end

        % Find Pz electrode 
        posteriorElectrodesIdx = [];
        posteriorElectrodesLabels = {};
        
        if isfield(S,'chanlocs') && isstruct(S.chanlocs) && isfield(S.chanlocs,'labels')
        
            labels = strtrim(string({S.chanlocs.labels}));
            
            datasetName = regexprep(dataMainFolders{iiDataset}, '_cond\d+_stim\d+_$', '');

            switch datasetName 
                case 'childStoriesEEG_Sysoeva'
                    posteriorElectrodesLabels = {'Pz','P3','P4','CP1','CP2'};
                case 'LalorNatSpeech'
                    posteriorElectrodesLabels = {'A3','A4','A19','A20','A31','A32'};
                case 'hyperscanningAmir_Listening'
                    posteriorElectrodesLabels = {'CPz','Pz','P3','P4', 'POz'};
                otherwise 
                    posteriorElectrodesLabels = {'Pz','P1','P2', 'CPz','CP1','CP2'};

        
            end

            posteriorElectrodesIdx = find(ismember(labels, string(posteriorElectrodesLabels)));
        end
        
        % Select N400 latency window
        samples = t >= n400Start & t <= n400End;
        
        if isempty(posteriorElectrodesIdx) || any(posteriorElectrodesIdx > size(w,2))
        
            warning('Posterior ROI unavailable for %s',filename);
        elseif ~any(samples)
            warning('N400 latency window unavailable: %s',filename);
        else
        
            %% Measuring the mean / peak negative amplitude
        
            posteriorTRF = mean(w(:, posteriorElectrodesIdx), 2, 'omitnan'); % avg TRF across selcted posterior electrodes
            n400ValRaw(iiDataset, ii) = mean(posteriorTRF(samples), 'omitnan'); % mean signed amplidude in the n400 window
            n400PeakRaw(iiDataset, ii) = min(posteriorTRF(samples));
        
        end
        %% Prediction correlation (full scalp)
        if ~isempty(rAllElec)
            pcorr_nonnorm(iiDataset,ii) = mean(rAllElec(:),'omitnan');

            validCorr = rAllElec(isfinite(rAllElec));

            if ~isempty(validCorr)
                p95 = prctile(validCorr,95);

                if p95 > 0
                    predcorr(iiDataset,ii) = pcorr_nonnorm(iiDataset,ii)/p95;
                end
            end
        end

        %% Prediction correlation at posterior electrodes
        if ~isempty(posteriorElectrodesIdx) && ~isempty(rAllElec)
        
            % rAllElec is expected to be channels x subjects
            if size(rAllElec,1) == numel(S.chanlocs)
        
                posteriorR = rAllElec(posteriorElectrodesIdx,:);
        
            elseif size(rAllElec,2) == numel(S.chanlocs)
        
                posteriorR = rAllElec(:,posteriorElectrodesIdx);
        
            else
        
                warning('Unexpected correlation dimensions: %s',filename);
                posteriorR = [];
        
            end
        
            if ~isempty(posteriorR)
                pcorr_posterior(iiDataset,ii) = mean(posteriorR(:),'omitnan');
            end
        
        end

        % Store number of processed chunks
        dataLen(iiDataset) = ii;

    end

    if dataLen(iiDataset) == 0
        disp(['Empty dataset: ',dataMainFoldersLabels{iiDataset}])

    else

        % Normalise N400 across time within each dataset
        valid = isfinite(n400ValRaw(iiDataset,:));

        if sum(valid) >= 2 && std(n400ValRaw(iiDataset,valid)) > 0

            n400Val(iiDataset,valid) = zscore(n400ValRaw(iiDataset,valid));
        end
        if sum(valid) >= 2 && std(n400PeakRaw(iiDataset,valid)) > 0

            n400Peak(iiDataset,valid) = zscore(n400PeakRaw(iiDataset,valid));
        end

        time = ((1:dataLen(iiDataset))-0.5)*defaultChunkLength;

        % N400 over time
        figure(10);

        subplot(1,2,1)
        hold on;
        plot(time, n400ValRaw(iiDataset,1:dataLen(iiDataset)), 'DisplayName',dataMainFoldersLabels{iiDataset});

        xlabel('Time (s)')
        ylabel('Mean signed Pz TRF amplitude')
        title('N400 window (300-500 ms) - Raw')

        subplot(1,2,2)
        hold on;
        plot(time, ...
            n400Val(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});

        xlabel('Time (s)')
        ylabel('Normalised signed amplitude')
        title('N400 window (300-500 ms) - Normalised')

        % N400 negative peak over time
        figure(12);
        
        subplot(1,2,1)
        hold on;
        plot(time, n400PeakRaw(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});
        
        xlabel('Time (s)')
        ylabel('Most negative Pz TRF weight')
        title('N400 window peak (300-500 ms) - Raw')
        
        subplot(1,2,2)
        hold on;
        plot(time, n400Peak(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});
        
        xlabel('Time (s)')
        ylabel('Normalised negative peak')
        title('N400 window peak (300-500 ms) - Normalised')

        % Prediction correlations over time
        figure(11);

        subplot(1,2,1)
        hold on;
        plot(time, pcorr_nonnorm(iiDataset,1:dataLen(iiDataset)), 'DisplayName',dataMainFoldersLabels{iiDataset});

        xlabel('Time (s)')
        ylabel('Mean prediction correlation')
        title('Raw prediction correlation')

        subplot(1,2,2)
        hold on;
        plot(time, predcorr(iiDataset,1:dataLen(iiDataset)), 'DisplayName',dataMainFoldersLabels{iiDataset});

        xlabel('Time (s)')
        ylabel('Normalised prediction correlation')
        title('Normalised prediction correlation')

    end
end

% Add legends after all datasets have been plotted
figure(10)
subplot(1,2,1)
legend('show','Location','best')
subplot(1,2,2)
legend('show','Location','best')

figure(11)
subplot(1,2,1)
legend('show','Location','best')
subplot(1,2,2)
legend('show','Location','best')

figure(12)
subplot(1,2,1)
legend('show','Location','best')
subplot(1,2,2)
legend('show','Location','best')

% Save calculated results
save(fullfile(resultsDir,'wordOnset_timeAnalysis_summary.mat'), ...
    'dataMainFolders','dataMainFoldersLabels', ...
    'n400Val','n400ValRaw', ...
    'predcorr','pcorr_nonnorm','dataLen', ...
    'defaultChunkLength','n400Start','n400End');

fprintf('\nAnalysis completed.\n');

for iiDataset = 1:nDatasets
    fprintf('%s: %d chunks\n', dataMainFoldersLabels{iiDataset}, dataLen(iiDataset));
end

%% Average N400 across datasets

% Use the first 12 chunks (0-120 seconds)
if superSub
    nPlotChunks = max(dataLen);
else
    nPlotChunks = 9; 
end

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Keep datasets with all 12 chunks and N400 measurements
validDatasets = dataLen >= nPlotChunks & all(isfinite(n400Val(:,1:nPlotChunks)),2);

n400Data = n400Val(validDatasets,1:nPlotChunks);

if isempty(n400Data)
    warning('No datasets with 12 valid N400 chunks.');
else

    % Average normalized signed amplitude across datasets
    meanN400 = mean(n400Data,1,'omitnan');
    stdN400 = std(n400Data, 0,1, 'omitnan');

    nN400 = sum(isfinite(n400Data), 1);
    semN400 = stdN400 ./ sqrt(nN400);

    figure(20);
    clf;
    hold on;
    
    plot(time,meanN400,'-o','LineWidth',2);
    errorbar(time, meanN400, semN400, '-o', 'LineWidth', 2);

    yline(0,'--');

    xlabel('Time (s)');
    ylabel('Normalised signed Pz amplitude');
    title('Average N400-window response across datasets');

    grid on;

    fprintf('Datasets included in N400 average: %d\n', sum(validDatasets));

end


% Keep datasets with at least 9 valid peak measurements
validPeakDatasets = dataLen >= nPlotChunks & all(isfinite(n400Peak(:,1:nPlotChunks)),2);

% Average N400 negative peak with error bars
peakData = n400Peak(validPeakDatasets,1:nPlotChunks);

if ~isempty(peakData)

    meanPeak = mean(peakData,1,'omitnan');
    stdPeak = std(peakData, 0, 1, 'omitnan');

    nPeak = sum(isfinite(peakData), 1);
    semPeak = stdPeak ./ sqrt(nPeak);

    figure(21);
    plot(time,meanPeak,'-o','LineWidth',2);
    errorbar(time, meanPeak, semPeak, '-o', 'LineWidth', 2);
    yline(0,'--');

    xlabel('Time (s)');
    ylabel('Normalised negative peak at the considered electrodes');
    title('Average N400-window negative peak across datasets');

    grid on;

    fprintf('Datasets included in peak average: %d\n', sum(validPeakDatasets));

end

%% Average prediction correlation across datasets

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Keep datasets with at least 9 valid chunks
validCorrDatasets = dataLen >= nPlotChunks & all(isfinite(pcorr_nonnorm(:,1:nPlotChunks)),2);

% Extract raw prediction correlations
corrData = pcorr_nonnorm(validCorrDatasets,1:nPlotChunks);

if isempty(corrData)

    warning('No datasets with 9 valid prediction correlation chunks.');

else

    % Average raw prediction correlations across datasets
    meanCorr = mean(corrData,1,'omitnan');
    stdCorr = std(corrData, 0, 1, 'omitnan');
    nCorr = sum(isfinite(corrData), 1);
    semCorr = stdCorr ./ sqrt(nCorr);

    figure(22);
    clf;
    hold on;
    
    errorbar(time, meanCorr, semCorr, '-o', 'LineWidth', 2);
    plot(time,meanCorr,'-o','LineWidth',2);

    xlabel('Time (s)');
    ylabel('Mean prediction correlation (r)');
    title('Average word-onset prediction correlation across datasets');

    grid on;

    fprintf('Datasets included in correlation average: %d\n', sum(validCorrDatasets));

end


%% Average prediction correlation across datasets (posterior)

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Keep datasets with at least 9 valid chunks
validCorrDatasets = dataLen >= nPlotChunks & all(isfinite(pcorr_posterior(:,1:nPlotChunks)),2);

% Extract raw prediction correlations
corrDataPosterior = pcorr_posterior(validCorrDatasets,1:nPlotChunks);

if isempty(corrDataPosterior)

    warning('No datasets with 9 valid prediction correlation chunks.');

else

    % Average raw prediction correlations across datasets
    meanCorrPos = mean(corrDataPosterior,1,'omitnan');
    stdCorrPos = std(corrDataPosterior, 0, 1, 'omitnan');
    nCorrPos = sum(isfinite(corrDataPosterior), 1);
    semCorrPos = stdCorrPos ./ sqrt(nCorrPos);

    figure(24);
    clf;
    hold on;
    
    errorbar(time, meanCorrPos, semCorrPos, '-o', 'LineWidth', 2);
    plot(time,meanCorrPos,'-o','LineWidth',2);

    xlabel('Time (s)');
    ylabel('Mean prediction correlation (r)');
    title('Average word-onset prediction correlation (of posterior electrodes) across datasets');

    grid on;

    fprintf('Datasets included in correlation average: %d\n', sum(validCorrDatasets));

end


% average signal 
figure;
plot(t, w(:,posteriorElectrodesIdx), 'LineWidth', 1.5);
xline(300,'--');
xline(500,'--');
yline(0,'--');

xlabel('Latency after word onset (ms)');
ylabel('TRF weight');
title('Word-onset TRF at the considered electrodes');
grid on;


%% Average word count across datasets

% Each dataset contributes its mean word count for each available chunk
wordCountData = nan(nDatasets, nChunks);

for iiDataset = 1:nDatasets

    for ii = 1:nChunks

        filename = "trfOut_100sChunk_" + ...
            dataMainFolders{iiDataset} + ii + ".mat";

        filepath = fullfile(resultsDir, filename);

        if ~isfile(filepath)
            continue
        end

        S = load(filepath, 'wordCountMeanAll');

        if ~isfield(S, 'wordCountMeanAll') || ...
                isempty(S.wordCountMeanAll)
            continue
        end

        % Average across subjects within this dataset
        wordCountData(iiDataset,ii) = ...
            mean(S.wordCountMeanAll, 'omitnan');

    end
end

%% Plot average words per chunk

% Include all datasets available for each chunk
meanWords = mean(wordCountData, 1, 'omitnan');
stdWords = std(wordCountData, 0, 1, 'omitnan');

% Number of datasets contributing to each chunk
nAvailable = sum(isfinite(wordCountData), 1);

% Time at the centre of each chunk
time = (0:nChunks-1)*defaultChunkLength + defaultChunkLength/2;

figure(23);
clf;

bar(time, meanWords);
hold on;

errorbar(time, meanWords, stdWords, ...
    'k.', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Average number of words');
title('Average word count per 10-second chunk');

xticks(time);
grid on;

% Show number of contributing datasets
fprintf('\nDatasets contributing to each word-count chunk:\n');
disp(nAvailable);

% Save results
save(fullfile(resultsDir, 'wordCount_summary.mat'), ...
    'wordCountData', 'meanWords', 'stdWords', ...
    'nAvailable', 'time');

%% Specific Chunks Analysis 
chunkOfInterest = 4;   % centred at 35 s

disp('Datasets contributing to the 35 s N400 average:')
disp(dataMainFoldersLabels(validDatasets))

disp('Normalised N400 values at 35 s:')
disp(n400Val(validDatasets, chunkOfInterest))

disp('Normalised peak values at 35 s:')
disp(n400Peak(validPeakDatasets, chunkOfInterest))


%% Plot individual dataset N400 traces

% Recreate the correct time vector for 9 chunks
timeN400 = (0:nPlotChunks-1)*defaultChunkLength ...
    + defaultChunkLength/2;

figure(30);
clf;
hold on;

for i = find(validDatasets)'

    plot(timeN400, n400Val(i,1:nPlotChunks), '-o', ...
        'DisplayName', dataMainFoldersLabels{i});

end

xline(35,'r--','LineWidth',1.5);
yline(0,'k--');

xlabel('Time (s)');
ylabel('Normalised signed Pz amplitude');
title('Individual dataset N400 traces');

legend('show','Location','bestoutside');
grid on;


%% Compare chunk 4 against neighbouring chunks

chunkOfInterest = 4;

datasetNames = string(dataMainFoldersLabels(validDatasets))';

values25 = n400Data(:,3);
values35 = n400Data(:,4);
values45 = n400Data(:,5);

% Difference from the average of neighbouring chunks
difference35 = values35 - (values25 + values45)/2;

results35 = table( ...
    datasetNames, values25, values35, values45, difference35, ...
    'VariableNames', ...
    {'Dataset','At25s','At35s','At45s','Difference35'});

% Sort from strongest positive change to strongest negative
results35 = sortrows(results35,'Difference35','descend');

disp(results35);


%% Leave-one-dataset-family-out N400 analysis

% Use the same 9 chunks and eligible entries as Figure 20
timeN400 = (0:nPlotChunks-1)*defaultChunkLength ...
    + defaultChunkLength/2;

names = string(dataMainFoldersLabels(validDatasets));
data = n400Val(validDatasets,1:nPlotChunks);

% Group conditions from the same dataset family
families = regexprep(names, '_cond\d+_stim\d+_$', '');

uniqueFamilies = unique(families, 'stable');

figure(34);
clf;
hold on;

% Original mean
plot(timeN400, mean(data,1,'omitnan'), ...
    'k-o', 'LineWidth', 2.5, ...
    'DisplayName', 'All datasets');

% Recompute mean after excluding each family
for k = 1:numel(uniqueFamilies)

    keep = families ~= uniqueFamilies(k);

    if ~any(keep)
        continue
    end

    meanExcluded = mean(data(keep,:),1,'omitnan');

    plot(timeN400, meanExcluded, '-o', ...
        'DisplayName', ['Without ' char(uniqueFamilies(k))]);

end

xline(35,'r--','HandleVisibility','off');
yline(0,'k--','HandleVisibility','off');

xlabel('Time (s)');
ylabel('Normalised signed Pz amplitude');
title('N400 sensitivity to dataset-family exclusion');

legend('show','Location','bestoutside');
grid on;


%% Average N400 giving equal weight to each dataset family

names = string(dataMainFoldersLabels(validDatasets));
data = n400Val(validDatasets,1:nPlotChunks);

families = regexprep(names, '_cond\d+_stim\d+_$', '');
uniqueFamilies = unique(families, 'stable');

familyData = nan(numel(uniqueFamilies),nPlotChunks);

for k = 1:numel(uniqueFamilies)
    idx = families == uniqueFamilies(k);
    familyData(k,:) = mean(data(idx,:),1,'omitnan');
end

meanFamily = mean(familyData,1,'omitnan');
semFamily = std(familyData,0,1,'omitnan') ...
    ./ sqrt(sum(isfinite(familyData),1));

timeN400 = (0:nPlotChunks-1)*defaultChunkLength ...
    + defaultChunkLength/2;

figure(35);
clf;
hold on;

errorbar(timeN400,meanFamily,semFamily,'-o','LineWidth',2);
yline(0,'--');
xline(35,'r--');

xlabel('Time (s)');
ylabel('Normalised signed Pz amplitude');
title('Average N400 across dataset families (equal weighting)');
grid on;

fprintf('Number of dataset families: %d\n',numel(uniqueFamilies));


%% Compare word count and N400 across dataset families

% Use the same entries included in the N400 analysis
wordData = wordCountData(validDatasets,1:nPlotChunks);

% Group word counts by dataset family
familyWordData = nan(numel(uniqueFamilies),nPlotChunks);

for k = 1:numel(uniqueFamilies)

    idx = families == uniqueFamilies(k);

    % Mean across conditions within each family
    familyWordData(k,:) = mean(wordData(idx,:),1,'omitnan');

end

% Average across dataset families
meanFamilyWords = mean(familyWordData,1,'omitnan');

nFamilyWords = sum(isfinite(familyWordData),1);

semFamilyWords = std(familyWordData,0,1,'omitnan') ...
    ./ sqrt(nFamilyWords);

% N400 values from the equal-family-weighting analysis
meanFamilyN400 = mean(familyData,1,'omitnan');

semFamilyN400 = std(familyData,0,1,'omitnan') ...
    ./ sqrt(sum(isfinite(familyData),1));

timeN400 = (0:nPlotChunks-1)*defaultChunkLength ...
    + defaultChunkLength/2;

figure(36);
clf;

% Word count
subplot(2,1,1);
bar(timeN400,meanFamilyWords);
hold on;

errorbar(timeN400,meanFamilyWords,semFamilyWords, ...
    'k.','LineWidth',1.5);

xline(35,'r--');
xlabel('Time (s)');
ylabel('Mean word count');
title('Average words per chunk (equal family weighting)');
xticks(timeN400);
grid on;

% N400
subplot(2,1,2);

errorbar(timeN400,meanFamilyN400,semFamilyN400, ...
    '-o','LineWidth',2);

xline(35,'r--');
yline(0,'k--');

xlabel('Time (s)');
ylabel('Normalised signed Pz amplitude');
title('Average N400 response (equal family weighting)');
xticks(timeN400);
grid on;


%% Compare word counts around 35 seconds

wordResults35 = table( ...
    uniqueFamilies(:), ...
    familyWordData(:,3), ...
    familyWordData(:,4), ...
    familyWordData(:,5), ...
    familyWordData(:,4) - ...
        (familyWordData(:,3)+familyWordData(:,5))/2, ...
    'VariableNames', ...
    {'Dataset','Words25s','Words35s','Words45s','Difference35'});

disp(wordResults35);

%% 
