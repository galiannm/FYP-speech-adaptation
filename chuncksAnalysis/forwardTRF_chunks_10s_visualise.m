
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

% dataMainFolders = {'LalorNatSpeech_cond1_stim1_/','LalorNatSpeechReverse_cond1_stim1_/','cocktailAttSwitch_cond1_stim1_/','cocktailAttSwitch_cond1_stim2_/','laughterEEG_cond1_stim1_/','diliBach_cond1_stim1_/','polyphonicBach_cond1_stim1_/','polyphonicBach_cond1_stim2_/','spectrotemporalModulationsIgnacio_cond1_stim1_/','spectrotemporalModulationsIgnacio_cond2_stim1_/','spectrotemporalModulationsIgnacio_cond3_stim1_/','spectrotemporalModulationsIgnacio_cond4_stim1_/','spectrotemporalModulationsIgnacio_cond5_stim1_/'};
% dataMainFoldersLabels = {'LalorNatSpeech','LalorNatSpeechReverse','cocktailAttSwitch-Target','cocktailAttSwitch-Masker','laughterEEG','diliBach','polyphonicBach-main','polyphonicBach-support','spectrotempMod-clean','spectrotempMod-temp60','spectrotempMod-temp80','spectrotempMod-spectr60','spectrotempMod-spectr80'};

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

% dataMainFolders = {'LalorNatSpeech/','LalorNatSpeechReverse/','diliBach/','laughterEEG/','cocktailAttSwitch/','cocktailAttSwitch/','spectrotemporalModulationsIgnacio/','polyphonicBach/','polyphonicBach/'};
% dataMainFoldersLabels = {'LalorNatSpeech','LalorNatSpeechReverse','diliBach','laughterEEG','cocktailAttSwitch_target','cocktailAttSwitch_masker','spectrotempMod-clean','spectrotempMod-temp60','spectrotempMod-temp80','spectrotempMod-spectr60','spectrotempMod-spectr80','polyphonicBachMain/','polyphonicBachSupport/'};
% % dataMainFolders = {'spectrotemporalModulationsIgnacio/'};
% % toPreprocess = [0];
% % dataMainFolders = {'musicImagery/'};
% toPreprocess = [0,0,0,0,0,0,0,1,0];
% stimFeatureDataset = [1,1,1,1,1,2,1,1,2];

% dataMainFolders = {'emilyDialogue/', 'anxietyEEG/'};
nChunks = 20;
clear n1Idx n1Val p2Idx p2Val
clear dataLen

% Initialise results for missing chunks
n1Val = nan(length(dataMainFolders), nChunks);
p2Val = nan(length(dataMainFolders), nChunks);
n1Idx = nan(length(dataMainFolders), nChunks);
p2Idx = nan(length(dataMainFolders), nChunks);
dataLen = zeros(length(dataMainFolders), 1);

% 'AliceSpeech/',
for iiDataset = 1:length(dataMainFolders)
    figure(1);
    cc = 1+(iiDataset-1)*nChunks;

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks

        filename = "trfOut_10sChunk_" + ...
            dataMainFolders{iiDataset} + ii + ".mat";

        filepath = fullfile(resultsDir, filename);

        if ~isfile(filepath)
            dataLen(iiDataset) = ii-1;
            fprintf('Missing chunk: %s\n', filename);
            break
        end

        load(filepath, ...
            'chanlocs','modelAll','rAllElec', ...
            'bandpassFilterRange','downFs','reRefType', ...
            'dataMainFolder','nSubs','stimIdx', ...
            'tmin','tmax','lambdas');

        if ii == nChunks
            dataLen(iiDataset) = ii;
        end

        % Plot average TRF
        plotNormFlag = 0;
        avgModel = mTRFmodelAvg(modelAll,plotNormFlag);
        avgModel.w = squeeze(avgModel.w);

        t = avgModel.t;

        % find latency N1
        tSt = 60;
        tEnd = 100;
        samples = t>=tSt & t<=tEnd;
        gfp = std(avgModel.w');
        gfp(~samples) = 0;
        [n1Val(iiDataset,ii),n1Idx(iiDataset,ii)] = max(gfp);

        % find latency P2
        tSt = 120;
        tEnd = 180;
        samples = t>=tSt & t<=tEnd;
        gfp = std(avgModel.w');
        gfp(~samples) = 0;
        [p2Val(iiDataset,ii),p2Idx(iiDataset,ii)] = max(gfp);

%         figure(2); % topos pred corr
        % subplot(length(dataMainFolders),nChunks,cc)
        % topoplot(mean(rAllElec,2),chanlocs,'electrodes','off')
        % caxis([-cmaxPred,cmaxPred])
        % title(num2str(ii))
        % cc = cc + 1;
        %         figure(2); % topos pred corr
        fprintf('Plotting %s, chunk %d: %d channels, chanlocs type = %s\n', ...
        dataMainFoldersLabels{iiDataset}, ii, ...
        size(rAllElec,1), class(chanlocs));
        subplot(length(dataMainFolders),nChunks,cc)

        predCorr = mean(rAllElec, 2);

        if isstruct(chanlocs) && ~isempty(chanlocs) && ...
                numel(chanlocs) == numel(predCorr)

            try
                topoplot(predCorr,chanlocs,'electrodes','off')
                caxis([-cmaxPred,cmaxPred])
            catch ME
                cla
                text(0.5,0.5,'Invalid chanlocs', ...
                    'HorizontalAlignment','center')
                axis off
                warning('Topoplot failed for %s chunk %d: %s', ...
                    dataMainFoldersLabels{iiDataset},ii,ME.message);
            end

        else
            % Cannot plot scalp topology without valid electrode locations
            text(0.5,0.5,'No chanlocs', ...
                'HorizontalAlignment','center')
            axis off
        end

        title(num2str(ii))
        cc = cc + 1;
    end

    if dataLen(iiDataset) == 0
        disp('Empty dataset')
    else
        n1Val(iiDataset,1:dataLen(iiDataset)) = ...
            zscore(n1Val(iiDataset,1:dataLen(iiDataset)));

        p2Val(iiDataset,1:dataLen(iiDataset)) = ...
            zscore(p2Val(iiDataset,1:dataLen(iiDataset)));

        figure(2);
        subplot(1,2,1)
        hold on; plot(n1Val(iiDataset,1:dataLen(iiDataset))')
        subplot(1,2,2)
        hold on; plot(n1Idx(iiDataset,1:dataLen(iiDataset))','.','MarkerSize',10)
        title('N1')

        figure(3);
        subplot(1,2,1)
        hold on; plot(p2Val(iiDataset,1:dataLen(iiDataset))')
        subplot(1,2,2)
        hold on; plot(p2Idx(iiDataset,1:dataLen(iiDataset))','.','MarkerSize',10)
        legend(dataMainFoldersLabels)
        title('P2')

        figure(10);
        subplot(1,2,1)
        % hold on;plot((1:dataLen(iiDataset))*defaultChunkLength,smooth(n1Val(iiDataset,1:dataLen(iiDataset))'))
        hold on;plot((1:dataLen(iiDataset))*defaultChunkLength,n1Val(iiDataset,1:dataLen(iiDataset))')
        % legend(dataMainFoldersLabels)
        xlabel('Time (s)')
        ylabel('Normalised magnitude')
        title('N1')

        subplot(1,2,2)
        % hold on;plot((1:dataLen(iiDataset))*defaultChunkLength,smooth(p2Val(iiDataset,1:dataLen(iiDataset))'))
        hold on;plot((1:dataLen(iiDataset))*defaultChunkLength,p2Val(iiDataset,1:dataLen(iiDataset))')
        legend(dataMainFoldersLabels)
        xlabel('Time (s)')
        ylabel('Normalised magnitude')
        title('P2')
    end

end

xlabel('Time chunk')
ylabel('Normalised P2 magnitude')

% Export figure if prepExport.m is available
if exist('prepExport.m', 'file')
    run prepExport.m
end

%%
xMean = nan(size(p2Val,2),1);
xStd = nan(size(p2Val,2),1);

for ii = 1:size(p2Val,2)
    tmp = n1Val(:,ii);
    tmp = tmp(isfinite(tmp));

    if ~isempty(tmp)
        xMean(ii) = mean(tmp);
        xStd(ii) = std(tmp)/sqrt(length(tmp));
    end
end

figure;
errorbar(((1:length(xMean))-1)*defaultChunkLength,xMean,xStd,'LineWidth',1.5)
xlim([-10,160])

hold on
xMean = nan(size(p2Val,2),1);
xStd = nan(size(p2Val,2),1);

for ii = 1:size(p2Val,2)
    tmp = p2Val(:,ii);
    tmp = tmp(isfinite(tmp));

    if ~isempty(tmp)
        xMean(ii) = mean(tmp);
        xStd(ii) = std(tmp)/sqrt(length(tmp));
    end
end

errorbar(((1:length(xMean))-1)*defaultChunkLength,xMean,xStd,'LineWidth',1.5)
legend({'N1','P2'})
xlabel('Time (s)')
ylabel('Normalised magnitude')

if exist('prepExport.m', 'file')
    run prepExport.m
end

%%
figure(4);
subplot(1,2,1)
bar(n1Idx)
subplot(1,2,2)
bar(p2Idx)

groupSize = 3;
aa1 = 1:groupSize;
aa2 = (groupSize+1):(2*groupSize);
figure(5);
n1IdxGrouped = [mean(n1Idx(:,aa1),2),mean(n1Idx(:,aa2),2)];
p2IdxGrouped = [mean(p2Idx(:,aa1),2),mean(p2Idx(:,aa2),2)];
subplot(1,2,1)
bar(n1IdxGrouped)
subplot(1,2,2)
bar(p2IdxGrouped)
