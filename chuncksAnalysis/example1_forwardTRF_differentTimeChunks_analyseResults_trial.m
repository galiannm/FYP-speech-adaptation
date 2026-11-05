%% Analyse results

clear all
close all

addpath ../libs/cnsp_utils
addpath ../libs/cnsp_utils/cnd
addpath ../libs/mTRF-Toolbox_v2/mtrf
addpath ../libs/NoiseTools
addpath ../libs/eeglab
eeglab

%% Analyse results

% defaultChunkLength = 20; % Seconds
close all
cmaxPred = 0.1;
dataMainFolders = {'cocktailAttSwitch_cond1_stim1_','cocktailAttSwitch_cond1_stim2_', 'diliBach_cond1_stim1_',...
    'LalorNatSpeech_cond1_stim1_','LalorNatSpeechReverse_cond1_stim1_','laughterEEG_cond1_stim1_',... 
    'spectrotemporalModulationsIgnacio_cond1_stim1_','spectrotemporalModulationsIgnacio_cond2_stim1_',...
    'spectrotemporalModulationsIgnacio_cond3_stim1_','spectrotemporalModulationsIgnacio_cond4_stim1_', ...
    'spectrotemporalModulationsIgnacio_cond5_stim1_', 'Speech-evoked_cond1_stim1_', 'Speech-evoked_cond2_stim1_'};
dataMainFoldersLabels = {'cocktailAttSwitch_cond1_stim1_','cocktailAttSwitch_cond1_stim2_', 'diliBach_cond1_stim1_',...
    'LalorNatSpeech_cond1_stim1_','LalorNatSpeechReverse_cond1_stim1_','laughterEEG_cond1_stim1_',... 
    'spectrotemporalModulationsIgnacio_cond1_stim1_','spectrotemporalModulationsIgnacio_cond2_stim1_',...
    'spectrotemporalModulationsIgnacio_cond3_stim1_','spectrotemporalModulationsIgnacio_cond4_stim1_', ...
    'spectrotemporalModulationsIgnacio_cond5_stim1_', 'Speech-evoked_cond1_stim1_', 'Speech-evoked_cond2_stim1_'};
% dataMainFolders = {'emilyDialogue/', 'anxietyEEG/'};
nChunks = 25; % max num of trials considered
clear n1Idx n1Val p2Idx p2Val
clear dataLen
% 'AliceSpeech/', 
for iiDataset = 1:length(dataMainFolders)
    figure(1);
    cc = 1+(iiDataset-1)*nChunks;

    % 0-10s, 10-20s, etc
    for ii = 1:nChunks
        try
            load("trialByTrialControl_results/trfOut_Trial_Chunk_"+dataMainFoldersLabels{iiDataset}+ii+".mat",'chanlocs','modelAll','rAllElec','bandpassFilterRange','downFs','reRefType','dataMainFolder','nSubs','stimIdx','tmin','tmax','lambdas')
            % load("resultsChunks/trfOut_20sChunk_"+dataMainFolders{iiDataset}(1:end-1)+ii+".mat",'chanlocs','modelAll','rAllElec','bandpassFilterRange','downFs','reRefType','dataMainFolder','nSubs','stimIdx','tmin','tmax','lambdas')
            if ii == nChunks, dataLen(iiDataset) = ii; end
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
        avgModel.w = squeeze(avgModel.w);
        
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
        
        
%         figure(2); % topos pred corr 
        subplot(length(dataMainFolders),nChunks,cc)
        topoplot(mean(rAllElec,2),chanlocs,'electrodes','off')
        caxis([-cmaxPred,cmaxPred])
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
