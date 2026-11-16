
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
defaultChunkLength = 10; % Seconds

% Exploratory posterior circle, matching the full-trial analysis.
posteriorRoiCentre = [0,-0.2];
posteriorRoiRadius = 0.15;
close all
cmaxPred = 0.1;

superSub = true;

% Path to saved results
scriptFolder = fileparts(mfilename('fullpath'));
addpath(scriptFolder); % Shared coordinate-based ROI selector

if superSub 
    resultsDir = fullfile(scriptFolder, 'resultsChunks', 'superSubject', 'wordOnset');
else
    resultsDir = fullfile(scriptFolder, 'resultsChunks', 'wordOnset');
end

% Dataset names matching the saved results
% Discover datasets automatically because word-onset feature indices
% can differ across datasets.
resultFiles = dir(fullfile(resultsDir, sprintf('trfOut_%dsChunk_*.mat', defaultChunkLength)));
dataMainFolders = {};

for iFile = 1:length(resultFiles)

    tokens = regexp(resultFiles(iFile).name, sprintf('^trfOut_%dsChunk_(.*)_(\\d+)\\.mat$', defaultChunkLength), 'tokens', 'once');

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

    fprintf('\nProcessing %s\n', dataMainFoldersLabels{iiDataset});

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks

        filename = "trfOut_" + defaultChunkLength + "sChunk_" + dataMainFolders{iiDataset} + ii + ".mat";

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

        % Use the same coordinate-based circle for every dataset.
        posteriorElectrodesIdx = selectScalpROI(S.chanlocs,size(w,2), ...
            posteriorRoiCentre,posteriorRoiRadius);
        if ii == 1
            fprintf('Posterior spatial ROI: %d electrodes\n',numel(posteriorElectrodesIdx));
            disp(string({S.chanlocs(posteriorElectrodesIdx).labels}));
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
        ylabel('Mean signed posterior ROI TRF amplitude')
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
        ylabel('Most negative posterior ROI TRF weight')
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
    'defaultChunkLength','n400Start','n400End', ...
    'posteriorRoiCentre','posteriorRoiRadius');

fprintf('\nAnalysis completed.\n');

for iiDataset = 1:nDatasets
    fprintf('%s: %d chunks\n', dataMainFoldersLabels{iiDataset}, dataLen(iiDataset));
end

%% Average N400 across datasets

% Plot the available recording span; contributors can differ at each chunk.
if superSub
    nPlotChunks = max(dataLen);
else
    nPlotChunks = 9; 
end

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Retain series with any valid chunk; omit missing entries at each time point.
validDatasets = any(isfinite(n400Val(:,1:nPlotChunks)),2);

n400Data = n400Val(validDatasets,1:nPlotChunks);
n400Data(~isfinite(n400Data)) = NaN;

if isempty(n400Data)
    warning('No valid N400 measurements available.');
else

    % Average normalized signed amplitude across datasets
    meanN400 = mean(n400Data,1,'omitnan');
    stdN400 = std(n400Data, 0,1, 'omitnan');

    nN400 = sum(isfinite(n400Data), 1);
    semN400 = stdN400 ./ sqrt(nN400);
    semN400(nN400 < 2) = NaN; % SEM needs at least two series.

    figure(20);
    clf;
    hold on;
    
    plot(time,meanN400,'-o','LineWidth',2);
    errorbar(time, meanN400, semN400, '-o', 'LineWidth', 2);

    yline(0,'--');

    xlabel('Time (s)');
    ylabel('Normalised signed posterior ROI amplitude');
    title('N400-window response: available series per chunk');

    grid on;

    fprintf('N400 contributors per chunk:\n'); disp(nN400);

end


% Include every series with at least one valid peak measurement
validPeakDatasets = any(isfinite(n400Peak(:,1:nPlotChunks)),2);

% Average N400 negative peak with error bars
peakData = n400Peak(validPeakDatasets,1:nPlotChunks);
peakData(~isfinite(peakData)) = NaN;

if ~isempty(peakData)

    meanPeak = mean(peakData,1,'omitnan');
    stdPeak = std(peakData, 0, 1, 'omitnan');

    nPeak = sum(isfinite(peakData), 1);
    semPeak = stdPeak ./ sqrt(nPeak);
    semPeak(nPeak < 2) = NaN; % SEM needs at least two series.

    figure(21);
    plot(time,meanPeak,'-o','LineWidth',2);
    errorbar(time, meanPeak, semPeak, '-o', 'LineWidth', 2);
    yline(0,'--');

    xlabel('Time (s)');
    ylabel('Normalised negative peak at the considered electrodes');
    title('N400-window negative peak: available series per chunk');

    grid on;

    fprintf('Peak contributors per chunk:\n'); disp(nPeak);

end

%% Average prediction correlation across datasets

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Include every series with at least one valid correlation
validCorrDatasets = any(isfinite(pcorr_nonnorm(:,1:nPlotChunks)),2);

% Extract raw prediction correlations
corrData = pcorr_nonnorm(validCorrDatasets,1:nPlotChunks);
corrData(~isfinite(corrData)) = NaN;

if isempty(corrData)

    warning('No valid prediction correlations available.');

else

    % Average raw prediction correlations across datasets
    meanCorr = mean(corrData,1,'omitnan');
    stdCorr = std(corrData, 0, 1, 'omitnan');
    nCorr = sum(isfinite(corrData), 1);
    semCorr = stdCorr ./ sqrt(nCorr);
    semCorr(nCorr < 2) = NaN; % SEM needs at least two series.

    figure(22);
    clf;
    hold on;
    
    errorbar(time, meanCorr, semCorr, '-o', 'LineWidth', 2);
    plot(time,meanCorr,'-o','LineWidth',2);

    xlabel('Time (s)');
    ylabel('Mean prediction correlation (r)');
    title('Word-onset prediction correlation: available series per chunk');

    grid on;

    fprintf('Full-scalp correlation contributors per chunk:\n'); disp(nCorr);

end


%% Average prediction correlation across datasets (posterior)

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Include every series with at least one valid correlation
validCorrDatasets = any(isfinite(pcorr_posterior(:,1:nPlotChunks)),2);

% Extract raw prediction correlations
corrDataPosterior = pcorr_posterior(validCorrDatasets,1:nPlotChunks);
corrDataPosterior(~isfinite(corrDataPosterior)) = NaN;

if isempty(corrDataPosterior)

    warning('No valid prediction correlations available.');

else

    % Average raw prediction correlations across datasets
    meanCorrPos = mean(corrDataPosterior,1,'omitnan');
    stdCorrPos = std(corrDataPosterior, 0, 1, 'omitnan');
    nCorrPos = sum(isfinite(corrDataPosterior), 1);
    semCorrPos = stdCorrPos ./ sqrt(nCorrPos);
    semCorrPos(nCorrPos < 2) = NaN; % SEM needs at least two series.

    figure(24);
    clf;
    hold on;
    
    errorbar(time, meanCorrPos, semCorrPos, '-o', 'LineWidth', 2);
    plot(time,meanCorrPos,'-o','LineWidth',2);

    xlabel('Time (s)');
    ylabel('Mean prediction correlation (r)');
    title('Posterior prediction correlation: available series per chunk');

    grid on;

    fprintf('Posterior correlation contributors per chunk:\n'); disp(nCorrPos);

end


%% Contributor counts for the four aggregate curves
% A series is a dataset/condition/feature entry, not an independent participant.
contributorCounts = [sum(isfinite(n400Data),1); ...
    sum(isfinite(peakData),1); sum(isfinite(corrData),1); ...
    sum(isfinite(corrDataPosterior),1)];
figure('Name','Contributors to each chunk average');
plot(time,contributorCounts','-o','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Number of contributing series');
legend('N400 mean','N400 peak','Full-scalp correlation','Posterior correlation', ...
    'Location','best');
title('Available contributors per chunk (composition changes over time)');
yticks(0:nDatasets);
grid on;
fprintf('Contributor counts (rows: mean, peak, full-scalp r, posterior r):\n');
disp(contributorCounts);
save(fullfile(resultsDir,'wordOnset_aggregate_summary.mat'), ...
    'time','contributorCounts','dataMainFoldersLabels', ...
    'n400Data','peakData','corrData','corrDataPosterior', ...
    'validDatasets','validPeakDatasets','posteriorRoiCentre','posteriorRoiRadius');

%% Average word count across datasets

% Each dataset contributes its mean word count for each available chunk
wordCountData = nan(nDatasets, nChunks);

for iiDataset = 1:nDatasets

    for ii = 1:nChunks

        filename = "trfOut_" + defaultChunkLength + "sChunk_" + dataMainFolders{iiDataset} + ii + ".mat";
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
disp(dataMainFoldersLabels(isfinite(n400Val(:,chunkOfInterest))));

disp('Normalised N400 values at 35 s:')
disp(n400Val(isfinite(n400Val(:,chunkOfInterest)), chunkOfInterest))

disp('Normalised peak values at 35 s:')
disp(n400Peak(isfinite(n400Peak(:,chunkOfInterest)), chunkOfInterest))


%% Plot individual dataset N400 traces

% Recreate the time vector for the chosen recording span.
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
ylabel('Normalised signed posterior ROI amplitude');
title('Individual dataset N400 traces');

legend('show','Location','bestoutside');
grid on;


%% Compare chunk 4 against neighbouring chunks

% chunkOfInterest = 4;
% 
% datasetNames = string(dataMainFoldersLabels(validDatasets))';
% 
% values25 = n400Data(:,3);
% values35 = n400Data(:,4);
% values45 = n400Data(:,5);
% 
% % Difference from the average of neighbouring chunks
% difference35 = values35 - (values25 + values45)/2;
% 
% results35 = table( ...
%     datasetNames, values25, values35, values45, difference35, ...
%     'VariableNames', ...
%     {'Dataset','At25s','At35s','At45s','Difference35'});
% 
% % Sort from strongest positive change to strongest negative
% results35 = sortrows(results35,'Difference35','descend');
% 
% disp(results35);


%% Leave-one-dataset-family-out N400 analysis

% Use the same span and available entries as Figure 20.
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
ylabel('Normalised signed posterior ROI amplitude');
title('N400 sensitivity to dataset-family exclusion');

legend('show','Location','bestoutside');
grid on;


%% Plot posterior TRF weights for each dataset
% Uses resultsDir, dataMainFolders and dataMainFoldersLabels
% already defined by this analysis script.

figure('Name','Posterior TRFs by dataset');
layout = tiledlayout('flow','TileSpacing','compact');

electrodeFig = figure('Name','Selected electrodes by dataset', ...
    'Units','normalized','Position',[0.05,0.08,0.90,0.80]);
electrodeLayout = tiledlayout(electrodeFig, ...
    ceil(numel(dataMainFolders)/2),2,'TileSpacing','compact','Padding','compact');
title(electrodeLayout,'Selected-electrode TRFs averaged across available chunks');

electrodeMapFig = figure('Name','Electrode locations matching TRF colours', ...
    'Units','normalized','Position',[0.05,0.08,0.90,0.80]);
electrodeMapLayout = tiledlayout(electrodeMapFig, ...
    ceil(numel(dataMainFolders)/2),2,'TileSpacing','compact','Padding','compact');
title(electrodeMapLayout,'Electrode colours match the corresponding TRF curves');

for d = 1:numel(dataMainFolders)
    curves = [];
    referenceT = [];
    electrodeSum = [];
    referenceElectrodeLabels = strings(0);
    referenceRoiIdx = [];
    referenceChanlocs = [];

    for chunk = 1:nChunks
        filename = sprintf('trfOut_10sChunk_%s%d.mat', ...
            dataMainFolders{d},chunk);
        filepath = fullfile(resultsDir,filename);

        if ~isfile(filepath)
            continue
        end

        S = load(filepath,'modelAll','chanlocs');
        % Average models without normalising their weights.
        model = mTRFmodelAvg(S.modelAll,0);
        lag = model.t(:);
        weights = reshape(model.w,numel(lag),[]);
        roiIdx = selectScalpROI(S.chanlocs,size(weights,2), posteriorRoiCentre,posteriorRoiRadius);

        electrodeLabels = strtrim(string({S.chanlocs.labels}));
        if isempty(referenceT)
            referenceT = lag;
            referenceElectrodeLabels = electrodeLabels;
            referenceRoiIdx = roiIdx;
            referenceChanlocs = S.chanlocs;
            electrodeSum = zeros(size(weights));
        elseif ~isequal(referenceT,lag) || ...
                ~isequal(referenceElectrodeLabels,electrodeLabels) || ...
                ~isequal(referenceRoiIdx,roiIdx)
            warning('Lag grid, channel order or ROI differs: skipping %s',filename);
            continue
        end
        electrodeSum = electrodeSum + weights;

        % One regional time course per chunk.
        curves(:,end+1) = mean(weights(:,roiIdx),2);
    end

    ax = nexttile(layout);
    hold(ax,'on');

    if isempty(curves)
        title(ax,dataMainFoldersLabels{d},'Interpreter','none');
        electrodeAx = nexttile(electrodeLayout,d);
        title(electrodeAx,dataMainFoldersLabels{d},'Interpreter','none');
        axis(electrodeAx,'off');
        mapAx = nexttile(electrodeMapLayout,d);
        title(mapAx,dataMainFoldersLabels{d},'Interpreter','none');
        axis(mapAx,'off');
        continue
    end

    % Grey: individual chunks. Blue: their average.
    plot(ax,referenceT,curves, ...
        'Color',[0.8,0.8,0.8],'HandleVisibility','off');
    plot(ax,referenceT,mean(curves,2), ...
        'Color',[0,0.447,0.741],'LineWidth',2);

    xline(ax,300,'k--','HandleVisibility','off');
    xline(ax,500,'k--','HandleVisibility','off');
    yline(ax,0,'k:','HandleVisibility','off');

    title(ax,sprintf('%s (%d chunks)', ...
        dataMainFoldersLabels{d},size(curves,2)), ...
        'Interpreter','none');
    xlabel(ax,'Latency after word onset (ms)');
    ylabel(ax,'Signed posterior TRF weight');
    grid(ax,'on');

    % One grid panel per series: a chunk-averaged curve per electrode.
    meanElectrodeWeights = electrodeSum / size(curves,2);
    electrodeAx = nexttile(electrodeLayout,d);
    electrodeLines = plot(electrodeAx,referenceT, ...
        meanElectrodeWeights(:,referenceRoiIdx),'LineWidth',1.5);
    electrodeColours = lines(numel(referenceRoiIdx));
    if numel(referenceRoiIdx) > 7
        electrodeColours = turbo(numel(referenceRoiIdx));
    end
    for e = 1:numel(electrodeLines)
        electrodeLines(e).Color = electrodeColours(e,:);
    end
    xline(electrodeAx,300,'--','HandleVisibility','off');
    xline(electrodeAx,500,'--','HandleVisibility','off');
    yline(electrodeAx,0,'--','HandleVisibility','off');
    xlabel(electrodeAx,'Latency after word onset (ms)');
    ylabel(electrodeAx,'Signed TRF weight');
    title(electrodeAx,sprintf('Selected electrodes: %s (mean of %d chunks)', ...
        dataMainFoldersLabels{d},size(curves,2)),'Interpreter','none');
    legend(electrodeAx,electrodeLines, ...
        referenceElectrodeLabels(referenceRoiIdx), ...
        'Interpreter','none','Location','bestoutside');
    grid(electrodeAx,'on');

    % Location-only topoplot: marker colour identifies an electrode, not weight.
    mapAx = nexttile(electrodeMapLayout,d);
    [~,locationXY] = selectScalpROI(referenceChanlocs, ...
        numel(referenceChanlocs),posteriorRoiCentre,posteriorRoiRadius);
    validLocations = all(isfinite(locationXY),2);
    locationRadii = sqrt(sum(locationXY(validLocations,:).^2,2));
    locationPlotRadius = max(0.5,min(1,max(locationRadii)*1.02));
    axes(mapAx);
    topoplot([],referenceChanlocs,'style','blank','electrodes','on', ...
        'plotrad',locationPlotRadius,'emarker',{'.',[0.65,0.65,0.65],5,1});
    hold(mapAx,'on');
    % topoplot plots [radius*sin(theta), radius*cos(theta)], scaled to 0.5.
    locationScale = 0.5/locationPlotRadius;
    locationMarkers = gobjects(numel(referenceRoiIdx),1);
    for e = 1:numel(referenceRoiIdx)
        ch = referenceRoiIdx(e);
        markerX = -locationXY(ch,1)*locationScale;
        markerY = locationXY(ch,2)*locationScale;
        locationMarkers(e) = plot3(mapAx,markerX,markerY,2.2,'o', ...
            'MarkerFaceColor',electrodeColours(e,:), ...
            'MarkerEdgeColor','k','MarkerSize',8);
    end
    legend(mapAx,locationMarkers,referenceElectrodeLabels(referenceRoiIdx), ...
        'Interpreter','none','Location','bestoutside');
    title(mapAx,dataMainFoldersLabels{d},'Interpreter','none');
end

%% Scalp topography for each dataset
figure('Name','300–500 ms TRF topographies by dataset');
layout = tiledlayout('flow','TileSpacing','compact');

for d = 1:numel(dataMainFolders)
    chunkMaps = [];
    referenceLabels = strings(0);
    mapChanlocs = [];

    for chunk = 1:nChunks
        filename = sprintf('trfOut_10sChunk_%s%d.mat', ...
            dataMainFolders{d},chunk);
        filepath = fullfile(resultsDir,filename);

        if ~isfile(filepath)
            continue
        end

        S = load(filepath,'modelAll','chanlocs');
        model = mTRFmodelAvg(S.modelAll,0);

        lag = model.t(:);
        weights = reshape(model.w,numel(lag),[]);
        window = lag >= 300 & lag <= 500;
        labels = strtrim(string({S.chanlocs.labels}));

        if ~any(window) || size(weights,2) ~= numel(labels)
            warning('Invalid dimensions or time window: %s',filename);
            continue
        end

        if isempty(chunkMaps)
            referenceLabels = labels;
            mapChanlocs = S.chanlocs;
        elseif ~isequal(referenceLabels,labels)
            warning('Channel order differs: skipping %s',filename);
            continue
        end

        % One 300–500 ms mean weight per electrode, per chunk.
        chunkMaps(:,end+1) = mean(weights(window,:),1)';
    end

    ax = nexttile(layout);

    if isempty(chunkMaps)
        title(ax,dataMainFoldersLabels{d},'Interpreter','none');
        axis(ax,'off');
        continue
    end

    % Average the available chunks with equal weight.
    mapValues = mean(chunkMaps,2);

    axes(ax);
    % Highlight exactly the circle-selected electrodes used in the curves.
    roiIdx = selectScalpROI(mapChanlocs,numel(mapValues), posteriorRoiCentre,posteriorRoiRadius);
    topoplot(mapValues,mapChanlocs, 'electrodes','on', 'emarker2',{roiIdx,'o','k',8,2});
    colorbar;

    title(ax,sprintf('%s (%d chunks)', ...
        dataMainFoldersLabels{d},size(chunkMaps,2)), ...
        'Interpreter','none');
end

title(layout,'Mean signed TRF weights: 300–500 ms');
