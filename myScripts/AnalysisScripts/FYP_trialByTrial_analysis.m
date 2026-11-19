scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));
addpath(scriptFolder); % Trial-violin plotting helper

addpath(fullfile(projectRoot,'libs','cnsp_utils'));
addpath(fullfile(projectRoot,'libs','cnsp_utils','cnd'));
addpath(fullfile(projectRoot,'libs','mTRF-Toolbox_v2','mtrf'));
addpath(fullfile(projectRoot,'libs','NoiseTools'));
addpath(fullfile(projectRoot,'libs','eeglab'));

eeglab

%% Analyse results

% defaultChunkLength = 20; % Seconds
close all
cmaxPred = 0.1;

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
resultsDir = fullfile(scriptFolder,'resultsTrials_envelope');

dataMainFoldersLabels = cell(size(dataMainFolders));

for d = 1:numel(dataMainFolders)
    % Automatically remove trailing folder separators.
    dataMainFolders{d} = regexprep(dataMainFolders{d},'[/\\]+$','');

    dataMainFoldersLabels{d} = sprintf('%s_stim%d', ...
        dataMainFolders{d},stimFeatureDataset(d));
end

dataMainFoldersLabels = cell(size(dataMainFolders));

for d = 1:numel(dataMainFolders)
    dataMainFoldersLabels{d} = sprintf('%s_stim%d', ...
        dataMainFolders{d},stimFeatureDataset(d));
end

nChunks = 25; % max num of trials considered
clear n1Idx n1Val p2Idx p2Val
clear dataLen

dataLen = zeros(1,numel(dataMainFolders));
n1Val = nan(numel(dataMainFolders),nChunks);
p2Val = nan(numel(dataMainFolders),nChunks);
n1Idx = nan(numel(dataMainFolders),nChunks);
p2Idx = nan(numel(dataMainFolders),nChunks);
predcorr = nan(numel(dataMainFolders),nChunks);

for iiDataset = 1:length(dataMainFolders)
    figure(1);
    cc = 1+(iiDataset-1)*nChunks;

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks
        try
            filename = sprintf('trfOut_Trial_%s_stim%d_trial%03d.mat', dataMainFolders{iiDataset},stimFeatureDataset(iiDataset),ii);
            filepath = fullfile(resultsDir,filename);

            load(filepath,'chanlocs','modelAll','rAllElec');
        catch
            dataLen(iiDataset) = ii-1;
            break
        end
        
        % Plot average TRF
        plotNormFlag = 0;
        isThereData = zeros(length(modelAll),1);
        for iiModel = 1:length(modelAll)
            isThereData(iiModel) = ~isempty(modelAll(iiModel).w);
        end
        modelAll = modelAll(find(isThereData));
        avgModel = mTRFmodelAvg(modelAll,plotNormFlag);

        if size(avgModel.w,1) ~= 1
            warning('%s: model has %d inputs; skipping single-envelope analysis.', ...
                dataMainFoldersLabels{iiDataset},size(avgModel.w,1));

            dataLen(iiDataset) = ii - 1;
            break
        end

        avgModel.w = reshape(avgModel.w,numel(avgModel.t),[]);
        
        t = avgModel.t;
        
        % find latency N1
        tSt = 50;
        tEnd = 120;
        samples = t>=tSt & t<=tEnd;
        gfp = std(avgModel.w');
        gfp(~samples) = 0;
        [n1Val(iiDataset,ii),n1Idx(iiDataset,ii)] = max(gfp);
            
        % find latency P2
        tSt = 120;
        tEnd = 220;
        samples = t>=tSt & t<=tEnd;
        gfp = std(avgModel.w');
        gfp(~samples) = 0;
        [p2Val(iiDataset,ii),p2Idx(iiDataset,ii)] = max(gfp);

        % find prediction correlations
       rAllElec(:,all(rAllElec == 0))=[]; % removes column if the entire column is zero
       predcorr(iiDataset, ii) = nanmean(rAllElec, "all")/prctile(rAllElec, 95, "all");
       dataLen(iiDataset) = ii;
        
        
        % figure(2); % topos pred corr 
        subplot(length(dataMainFolders),nChunks,cc)

        if isstruct(chanlocs) && ~isempty(chanlocs)
            % Convert coordinate storage types without changing their values.
            plotChanlocs = chanlocs;

            for ch = 1:numel(plotChanlocs)
                plotChanlocs(ch).theta = double(plotChanlocs(ch).theta);
                plotChanlocs(ch).radius = double(plotChanlocs(ch).radius);
            end

            topoR = mean(rAllElec,2,'omitnan');

            if any(isfinite(topoR))
                topoplot(topoR,plotChanlocs,'electrodes','off');
                caxis([-cmaxPred,cmaxPred]);
            else
                axis off
                text(0.5,0.5,'No valid prediction correlations', ...
                    'HorizontalAlignment','center');
            end
            caxis([-cmaxPred,cmaxPred]);
        else
            axis off
            text(0.5,0.5,'No electrode locations', ...
                'HorizontalAlignment','center');
        end

        title(num2str(ii))
        cc = cc + 1;
    end

    if dataLen(iiDataset) == 0
        disp('Empty dataset')
    else
        n1Val(iiDataset,1:dataLen(iiDataset)) = zscore(n1Val(iiDataset,1:dataLen(iiDataset)));
        p2Val(iiDataset,1:dataLen(iiDataset)) = zscore(p2Val(iiDataset,1:dataLen(iiDataset)));

%         figure(2); 
%         subplot(1,2,1)
%         hold on; plot(n1Val(iiDataset,1:dataLen(iiDataset))')
%         subplot(1,2,2)
%         hold on; plot(n1Idx(iiDataset,1:dataLen(iiDataset))')
% 
%         figure(3); 
%         subplot(1,2,1)
%         hold on; plot(p2Val(iiDataset,1:dataLen(iiDataset))')
%         subplot(1,2,2)
%         hold on; plot(p2Idx(iiDataset,1:dataLen(iiDataset))')
%         legend(dataMainFoldersLabels)
% 
%         figure(10); 
%         subplot(1,2,1)
%         hold on;plot((n1Val(iiDataset,1:dataLen(iiDataset))'))
%         % hold on;plot((n1Val(iiDataset,1:dataLen(iiDataset))'))
%         title('N1')
%         subplot(1,2,2)
%         hold on;plot((p2Val(iiDataset,1:dataLen(iiDataset))'))
%         % hold on;plot((p2Val(iiDataset,1:dataLen(iiDataset))'))
%         title('P2')
%         legend(dataMainFoldersLabels)
%         xlabel('Trial')
%         ylabel('Normalised magnitude')
     end
% 
 end
% 
% xlabel('Time chunk')
% ylabel('Normalised P2 magnitude')
% run prepExport.m
% 
% %%
% xMean = zeros(size(p2Val,2),1);
% xStd = zeros(size(p2Val,2),1);
% for ii = 1:size(p2Val,2)
%     tmp = n1Val(:,ii);
%     xMean(ii) = mean(tmp(tmp~=0));
%     xStd(ii) = std(tmp(tmp~=0))/sqrt(sum(tmp~=0));
% end
% 
% figure;
% errorbar(xMean,xStd,'LineWidth',1.5)
% xlim([-1,26])
% 
% hold on
% xMean = zeros(size(p2Val,2),1);
% xStd = zeros(size(p2Val,2),1);
% for ii = 1:size(p2Val,2)
%     tmp = p2Val(:,ii);
%     xMean(ii) = mean(tmp(tmp~=0));
%     xStd(ii) = std(tmp(tmp~=0))/sqrt(sum(tmp~=0));
% end
% errorbar(xMean,xStd,'LineWidth',1.5)
% legend({'N1','P2'})
% xlabel('Trial')
% ylabel('Normalised magnitude')
% run prepExport.m
% 
% %%
% figure(4);
% subplot(1,2,1)
% bar(n1Idx)
% subplot(1,2,2)
% bar(p2Idx)
% 
% 
% groupSize = 3;
% aa1 = 1:groupSize;
% aa2 = (groupSize+1):(2*groupSize);
% figure(5);
% n1IdxGrouped = [mean(n1Idx(:,aa1),2),mean(n1Idx(:,aa2),2)];
% p2IdxGrouped = [mean(p2Idx(:,aa1),2),mean(p2Idx(:,aa2),2)];
% subplot(1,2,1)
% bar(n1IdxGrouped)
% subplot(1,2,2)
% bar(p2IdxGrouped)

%% Violin summaries across dataset/feature series
% Keep the reference measurements and normalisation; change presentation only.
% Explicitly mask positions not processed for a series (including skipped series).
n1Violin = n1Val;
p2Violin = p2Val;
corrViolin = predcorr;
for d = 1:numel(dataMainFolders)
    n1Violin(d,dataLen(d)+1:end) = NaN;
    p2Violin(d,dataLen(d)+1:end) = NaN;
    corrViolin(d,dataLen(d)+1:end) = NaN;
end
[trialViolinFig,violinContributorCounts] = plotTrialViolins(n1Violin,p2Violin,corrViolin,dataMainFoldersLabels);
% No significance markers: inferential tests are a separate analysis.
