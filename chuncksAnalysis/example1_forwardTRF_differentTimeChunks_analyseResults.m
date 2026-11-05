
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
resultsDir = fullfile(scriptFolder, 'resultsChunks');

% Dataset names matching the saved results
dataMainFolders = { ...
    'LalorNatSpeech_cond1_stim1_', ...
    'LalorNatSpeechReverse_cond1_stim1_', ...
    'cocktailAttSwitch_cond1_stim1_', ...
    'cocktailAttSwitch_cond1_stim2_', ...
    'laughterEEG_cond1_stim1_', ...
    'diliBach_cond1_stim1_', ...
    'polyphonicBach_cond1_stim1_', ...
    'polyphonicBach_cond1_stim2_', ...
    'spectrotemporalModulationsIgnacio_cond1_stim1_', ...
    'spectrotemporalModulationsIgnacio_cond2_stim1_', ...
    'spectrotemporalModulationsIgnacio_cond3_stim1_', ...
    'spectrotemporalModulationsIgnacio_cond4_stim1_', ...
    'spectrotemporalModulationsIgnacio_cond5_stim1_'};

dataMainFoldersLabels = { ...
    'LalorNatSpeech', ...
    'LalorNatSpeechReverse', ...
    'cocktailAttSwitch-Target', ...
    'cocktailAttSwitch-Masker', ...
    'laughterEEG', ...
    'diliBach', ...
    'polyphonicBach-main', ...
    'polyphonicBach-support', ...
    'spectrotempMod-clean', ...
    'spectrotempMod-temp60', ...
    'spectrotempMod-temp80', ...
    'spectrotempMod-spectr60', ...
    'spectrotempMod-spectr80'};

nChunks = 20;
clear n1Idx n1Val p2Idx p2Val
clear dataLen

% Initialise results
nDatasets = length(dataMainFolders);

n1Val = nan(nDatasets,nChunks);
p2Val = nan(nDatasets,nChunks);
n1Idx = nan(nDatasets,nChunks);
p2Idx = nan(nDatasets,nChunks);

predcorr = nan(nDatasets,nChunks);
pcorr_nonnorm = nan(nDatasets,nChunks);

dataLen = zeros(nDatasets,1);

% 'AliceSpeech/',
for iiDataset = 1:length(dataMainFolders)

    cc = 1+(iiDataset-1)*nChunks;

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks

        filename = "trfOut_10sChunk_" + ...
            dataMainFolders{iiDataset} + ii + ".mat";

        filepath = fullfile(resultsDir,filename);

        if ~isfile(filepath)
            dataLen(iiDataset) = ii-1;
            fprintf('Missing chunk: %s\n',filename);
            break
        end

        load(filepath, ...
            'chanlocs','modelAll','rAllElec', ...
            'bandpassFilterRange','downFs','reRefType', ...
            'dataMainFolder','nSubs','stimIdx', ...
            'tmin','tmax','lambdas');

        % Plot average TRF
        plotNormFlag = 0;

        isThereData = zeros(length(modelAll),1);
        for iiModel = 1:length(modelAll)
            isThereData(iiModel) = ~isempty(modelAll(iiModel).w);
        end

        modelAll = modelAll(find(isThereData));

        if isempty(modelAll) || isempty(rAllElec)
            fprintf('Empty models: %s\n',filename);
            dataLen(iiDataset) = ii-1;
            break
        end

        avgModel = mTRFmodelAvg(modelAll,plotNormFlag);
        avgModel.w = squeeze(avgModel.w);

        t = avgModel.t;

        % find latency N1
        tSt = 60;
        tEnd = 100;
        samples = t>=tSt & t<=tEnd;
        gfp = std(avgModel.w');
        gfp(~samples) = 0;
        % [n1Val(iiDataset,ii),n1Idx(iiDataset,ii)] = max(gfp);
        n1Val(iiDataset,ii) = mean(gfp(samples));

        % find latency P2
        tSt = 120;
        tEnd = 180;
        samples = t>=tSt & t<=tEnd;
        gfp = std(avgModel.w');
        gfp(~samples) = 0;
        [p2Val(iiDataset,ii),p2Idx(iiDataset,ii)] = max(gfp);

        % find prediction correlations
        rAllElec(:,all(rAllElec == 0)) = []; % removes column if the entire column is zero

        if ~isempty(rAllElec)
            pcorr_nonnorm(iiDataset,ii) = mean(rAllElec(:),'omitnan');

            validCorr = rAllElec(isfinite(rAllElec));

            if ~isempty(validCorr)
                p95 = prctile(validCorr,95);

                if p95 > 0
                    predcorr(iiDataset,ii) = ...
                        pcorr_nonnorm(iiDataset,ii)/p95;
                end
            end
        end

        % Store number of processed chunks
        dataLen(iiDataset) = ii;

%         figure(2); % topos pred corr
%         subplot(length(dataMainFolders),nChunks,cc)
%         topoplot(mean(rAllElec,2),chanlocs,'electrodes','off')
%         caxis([-cmaxPred,cmaxPred])
%         title(num2str(ii))
%         cc = cc + 1;

    end

    if dataLen(iiDataset) == 0
        disp(['Empty dataset: ',dataMainFoldersLabels{iiDataset}])
    else

        n1Val(iiDataset,1:dataLen(iiDataset)) = ...
            zscore(n1Val(iiDataset,1:dataLen(iiDataset)));

        p2Val(iiDataset,1:dataLen(iiDataset)) = ...
            zscore(p2Val(iiDataset,1:dataLen(iiDataset)));

        % N1/P2 over time
        time = ((1:dataLen(iiDataset))-0.5)*defaultChunkLength;

        figure(10);

        subplot(1,2,1)
        hold on;
        plot(time, ...
            n1Val(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});
        xlabel('Time (s)')
        ylabel('Normalised magnitude')
        title('N1')

        subplot(1,2,2)
        hold on;
        plot(time, ...
            p2Val(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});
        xlabel('Time (s)')
        ylabel('Normalised magnitude')
        title('P2')

        % Prediction correlations over time
        figure(11);

        subplot(1,2,1)
        hold on;
        plot(time, ...
            pcorr_nonnorm(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});
        xlabel('Time (s)')
        ylabel('Mean prediction correlation')
        title('Raw prediction correlation')

        subplot(1,2,2)
        hold on;
        plot(time, ...
            predcorr(iiDataset,1:dataLen(iiDataset)), ...
            'DisplayName',dataMainFoldersLabels{iiDataset});
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

% Save calculated results
save(fullfile(resultsDir,'timeAnalysis_summary.mat'), ...
    'dataMainFolders','dataMainFoldersLabels', ...
    'n1Val','p2Val','n1Idx','p2Idx', ...
    'predcorr','pcorr_nonnorm','dataLen', ...
    'defaultChunkLength');

fprintf('\nAnalysis completed.\n');

for iiDataset = 1:nDatasets
    fprintf('%s: %d chunks\n', ...
        dataMainFoldersLabels{iiDataset}, ...
        dataLen(iiDataset));
end

%% Average N1/P2 across datasets
% tried to reproduce a similar figure as in the the original paper

% Use the first 12 chunks (0-120 seconds)
nPlotChunks = 12;

time = (0:nPlotChunks-1)*10 + 5;

% Keep datasets with all 12 chunks available
validDatasets = dataLen >= nPlotChunks;

n1Data = n1Val(validDatasets,1:nPlotChunks);
p2Data = p2Val(validDatasets,1:nPlotChunks);

% Average normalized magnitude across datasets
meanN1 = mean(n1Data,1,'omitnan');
meanP2 = mean(p2Data,1,'omitnan');

figure(20);
hold on;

plot(time,meanN1,'-o','LineWidth',2);
plot(time,meanP2,'-o','LineWidth',2);

xlabel('Time (s)');
ylabel('Normalised magnitude');
title('Average N1/P2 across datasets');

legend({'N1','P2'});
grid on;