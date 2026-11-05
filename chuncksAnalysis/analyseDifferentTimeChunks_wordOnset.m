
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
close all
cmaxPred = 0.1;

% Path to saved results
scriptFolder = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptFolder, 'resultsChunks', 'wordOnset');

% Dataset names matching the saved results
% Discover datasets automatically because word-onset feature indices
% can differ across datasets.
resultFiles = dir(fullfile(resultsDir, 'trfOut_10sChunk_*.mat'));

dataMainFolders = {};

for iFile = 1:length(resultFiles)

    tokens = regexp(resultFiles(iFile).name, ...
        '^trfOut_10sChunk_(.*)_(\d+)\.mat$', 'tokens', 'once');

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

dataLen = zeros(nDatasets,1);

% N400 latency window
n400Start = 300;
n400End = 500;

for iiDataset = 1:length(dataMainFolders)

    fprintf('\nProcessing %s\n', ...
        dataMainFoldersLabels{iiDataset});

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks

        filename = "trfOut_10sChunk_" + ...
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
        pzIdx = [];
        
        if isfield(S,'chanlocs') && isstruct(S.chanlocs) && isfield(S.chanlocs,'labels')
        
            labels = strtrim(string({S.chanlocs.labels}));
            pzIdx = find(strcmpi(labels,'Pz'),1);
        
        end
        
        % Select N400 latency window
        samples = t >= n400Start & t <= n400End;
        
        if isempty(pzIdx) || pzIdx > size(w,2)
        
            warning('Pz not found for %s',filename);
        
        else
        
            %% Measuring the mean / peak negative amplitude
        
            % Extract signed Pz TRF
            pzTRF = w(:,pzIdx);
        
            % Mean amplitude in the N400 latency window
            n400ValRaw(iiDataset,ii) = mean(pzTRF(samples),'omitnan');
            
            % Mean neg amplitude in the n400 window
            n400PeakRaw(iiDataset, ii) = min(pzTRF(samples));
        
        end

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
nPlotChunks = 9;

time = (0:nPlotChunks-1)*defaultChunkLength + defaultChunkLength/2;

% Keep datasets with all 12 chunks and N400 measurements
validDatasets = dataLen >= nPlotChunks & all(isfinite(n400Val(:,1:nPlotChunks)),2);

n400Data = n400Val(validDatasets,1:nPlotChunks);

if isempty(n400Data)

    warning('No datasets with 12 valid N400 chunks.');

else

    % Average normalized signed amplitude across datasets
    meanN400 = mean(n400Data,1,'omitnan');

    figure(20);
    hold on;

    plot(time,meanN400,'-o','LineWidth',2);

    yline(0,'--');

    xlabel('Time (s)');
    ylabel('Normalised signed Pz amplitude');
    title('Average N400-window response across datasets');

    grid on;

    fprintf('Datasets included in N400 average: %d\n', sum(validDatasets));

end


% Keep datasets with at least 9 valid peak measurements
validPeakDatasets = dataLen >= nPlotChunks & all(isfinite(n400Peak(:,1:nPlotChunks)),2);

peakData = n400Peak(validPeakDatasets,1:nPlotChunks);

if ~isempty(peakData)

    meanPeak = mean(peakData,1,'omitnan');

    figure(21);
    plot(time,meanPeak,'-o','LineWidth',2);
    yline(0,'--');

    xlabel('Time (s)');
    ylabel('Normalised negative peak at Pz');
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

    figure(22);
    hold on;

    plot(time,meanCorr,'-o','LineWidth',2);

    xlabel('Time (s)');
    ylabel('Mean prediction correlation (r)');
    title('Average word-onset prediction correlation across datasets');

    grid on;

    fprintf('Datasets included in correlation average: %d\n', sum(validCorrDatasets));

end


% average signal 
figure;
plot(t, w(:,pzIdx), 'LineWidth', 1.5);
xline(300,'--');
xline(500,'--');
yline(0,'--');

xlabel('Latency after word onset (ms)');
ylabel('TRF weight');
title('Word-onset TRF at Pz');
grid on;