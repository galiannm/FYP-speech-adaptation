
% Cognition and Natural Sensory Processing (CNSP) open-science initiative
% Forward TRF script
%
% Original script:
% CNSP2025_forwardTRF_example1.m
% Author: Giovanni M. Di Liberto (2025)
%
% FYP adaptation:
% Full-trial Word Onset TRF analysis (no 10-second chunking)

clear
close all
clc

%% ============ SETUP ============

% This script should be saved inside standardAnalysisScripts/
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = '/Users/manon/Documents/Uni/Trinity4/FYP/CNSP-resources/CNSP';

addpath(fullfile(projectRoot,'exampleScripts'));
addpath(fullfile(projectRoot,'standardAnalysisScripts'));
addpath(fullfile(projectRoot,'libs','cnsp_utils'));
addpath(fullfile(projectRoot,'libs','cnsp_utils','cnd'));
addpath(fullfile(projectRoot,'libs','mTRF-Toolbox_v2','mtrf'));
addpath(fullfile(projectRoot,'libs','NoiseTools'));
addpath(fullfile(projectRoot,'libs','eeglab'));

eeglab
close

%% ============ PARAMETERS ============

% Dataset and feature
datasetName = 'speechEmotionAdultEEG';
stimIdx = 3;              % Word Onset for LalorNatSpeech
condIdx = 1;

reRefType = 'Avg';
bandpassFilterRange = [1,8];

% TRF hyperparameters
dirTRF = 1;
tmin = -100;
tmax = 600;

lambdas = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4];

% Dataset locations
dataMainFolder = fullfile(projectRoot,'datasets2');
dataCNDfolder = fullfile(dataMainFolder,datasetName,'dataCND');
eegFolder = fullfile(dataCNDfolder,'1-8Hz');

% Results folder
resultDir = fullfile(projectRoot,'myOutputs', ...
    'fullTrial',datasetName,'wordOnset');

if ~exist(resultDir,'dir')
    mkdir(resultDir);
end

%% ============ FIND SUBJECTS ============

% Use already preprocessed EEG
eegFilenames = dir(fullfile(eegFolder,'pre_dataSub*.mat'));

% Exclude the previously generated super-subject
eegFilenames = eegFilenames( ...
    ~strcmp({eegFilenames.name},'pre_dataSub100.mat'));

nSubs = length(eegFilenames);

fprintf('Dataset: %s\n',datasetName);
fprintf('Found %d subjects\n',nSubs);

if nSubs == 0
    error('No preprocessed EEG files found in %s',eegFolder);
end

%% ============ LOAD STIMULUS ============

stimFilenameSame = fullfile(dataCNDfolder,'dataStim.mat');

if ~isfile(stimFilenameSame)
    error('Expected shared stimulus file: %s',stimFilenameSame);
end

load(stimFilenameSame,'stim');

fprintf('Loaded stimulus: %s\n',stimFilenameSame);

% Default to one condition if condition indices are absent
if ~isfield(stim,'condIdxs') || isempty(stim.condIdxs)
    stim.condIdxs = ones(1,size(stim.data,2));
end

stim.condIdxs = stim.condIdxs(:)';

fprintf('Feature selected: %d\n',stimIdx);
disp(stim.names(stimIdx));

%% ============ FORWARD TRF ============

rAll = [];
rAllElec = [];
modelAll = struct([]);
includedSubjectIdx = [];
bestLambdas = [];
chanlocs = [];

figure('Name','Full-trial TRF');

for sub = 1:nSubs

    %% Load preprocessed EEG
    eegPreFilename = fullfile(eegFolder,eegFilenames(sub).name);

    fprintf('\nProcessing %s\n',eegFilenames(sub).name);

    load(eegPreFilename,'eeg');

    subIdx = regexp(eegFilenames(sub).name,'\d+','match');
    subIdx = str2double(subIdx{1});

    %% Re-reference EEG
    eeg = cndReref(eeg,reRefType);

    %% Select stimulus feature
    stimFeature = stim;
    stimFeature.data = stimFeature.data(stimIdx,:);

    %% Select condition
    trials2keep = stimFeature.condIdxs == condIdx;

    if numel(eeg.data) ~= numel(stimFeature.data)
        warning('Subject %d: EEG and STIM trial counts differ',subIdx);
        continue
    end

    stimFeature.data = stimFeature.data(:,trials2keep);
    eeg.data = eeg.data(trials2keep);

    if numel(eeg.data) < 2
        warning('Subject %d: fewer than 2 trials',subIdx);
        continue
    end

    %% Match sampling frequencies
    % This preserves impulse timing for Word Onset features.

    eeg.fs = double(eeg.fs);
    stimFeature.fs = double(stimFeature.fs);

    if stimFeature.fs > eeg.fs
        fprintf('Downsampling STIM: %g -> %g Hz\n', ...
            stimFeature.fs,eeg.fs);

        stimFeature = cndDownsample(stimFeature,eeg.fs);
    end

    if stimFeature.fs ~= eeg.fs
        error('EEG and STIM sampling frequencies do not match');
    end

    %% Match EEG and stimulus lengths
    % Keep complete trials rather than fixed-length chunks.

    for tr = 1:numel(eeg.data)

        stimLen = size(stimFeature.data{tr},1);
        eegLen = size(eeg.data{tr},1);

        minLen = min(stimLen,eegLen);

        stimFeature.data{tr} = double( ...
            stimFeature.data{tr}(1:minLen,:));

        eeg.data{tr} = double(eeg.data{tr}(1:minLen,:));
    end

    %% TRF cross-validation
    % Use default ridge regularisation, matching CNSP2025_forwardTRF_example1.
    disp('Running mTRFcrossval');

    [stats,t] = mTRFcrossval( ...
        stimFeature.data, ...
        eeg.data, ...
        eeg.fs, ...
        dirTRF, ...
        tmin,tmax, ...
        lambdas, ...
        'verbose',0);

    % Select optimal regularisation parameter
    [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));

    fprintf('Prediction correlation r = %.5f\n',maxR);
    fprintf('Best lambda = %g\n',lambdas(bestLambda));

    %% Fit final TRF model
    % Keep Tikhonov for the final fit, as in CNSP2025_forwardTRF_example1.
    disp('Running mTRFtrain');

    model = mTRFtrain( ...
        stimFeature.data, ...
        eeg.data, ...
        eeg.fs, ...
        dirTRF, ...
        tmin,tmax, ...
        lambdas(bestLambda), ...
        'method','Tikhonov', ...
        'verbose',0);

    %% Store model
    if isempty(modelAll)
        modelAll = model;
        chanlocs = eeg.chanlocs;
    else
        modelAll(end+1) = model;
    end

    rAll(end+1) = maxR;

    rAllElec(:,end+1) = squeeze( ...
        mean(stats.r(:,bestLambda,:),1));

    includedSubjectIdx(end+1) = subIdx;
    bestLambdas(end+1) = lambdas(bestLambda);

    %% TRF visualisation
    mTRF_plotForwardTRF(eeg,modelAll,rAllElec);

    subplot(2,2,3);
    title('Word Onset TRF');

    drawnow;

    fprintf('Mean r so far = %.5f\n',mean(rAll));
end

%% ============ SAVE RESULTS ============

if isempty(modelAll)
    error('No TRF models were fitted.');
end

nIncluded = numel(modelAll);

save(fullfile(resultDir,'TRF_fullTrial.mat'), ...
    'modelAll','rAll','rAllElec', ...
    'chanlocs','includedSubjectIdx','bestLambdas', ...
    'datasetName','stimIdx','condIdx', ...
    'tmin','tmax','lambdas','bandpassFilterRange','reRefType');

savefig(gcf,fullfile(resultDir,'TRF_overview.fig'));

fprintf('\nSuccessfully fitted %d subjects\n',nIncluded);

%% ============ GROUP-LEVEL TRF ============

% Average the fitted TRF weights across subjects.
% No averaging of raw EEG across subjects.

t = modelAll(1).t;

nLags = numel(t);
nChans = size(modelAll(1).w,3);

W = zeros(nIncluded,nLags,nChans);

for s = 1:nIncluded
    W(s,:,:) = modelAll(s).w;
end

% Mean TRF across subjects
avgTRF = squeeze(mean(W,1));

%% ============ POSTERIOR SPATIAL ROI ============

% Circular ROI in the 2-D scalp projection (nose up, head radius = 0.5).
% x: left/right; y: front/back. Negative y selects the back of the head.
% These are geometric settings, independent of the fitted TRF amplitudes.
roiCentre = [0,-0.30];
roiRadius = 0.18;

if numel(chanlocs) ~= nChans || ...
        ~all(isfield(chanlocs,{'theta','radius'}))
    error('Spatial ROI requires theta/radius coordinates for every model channel.');
end
scalpXY = nan(nChans,2);
for ch = 1:nChans
    theta = chanlocs(ch).theta;
    radius = chanlocs(ch).radius;
    if isnumeric(theta) && isscalar(theta) && isfinite(theta) && ...
            isnumeric(radius) && isscalar(radius) && isfinite(radius) && radius >= 0
        scalpXY(ch,:) = [-radius*sind(theta),radius*cosd(theta)];
    end
end
validCoordinates = all(isfinite(scalpXY),2);
roiDistance = sqrt(sum((scalpXY-roiCentre).^2,2));
roiIdx = find(validCoordinates & roiDistance <= roiRadius);
if isempty(roiIdx)
    error('No electrodes fall within the spatial ROI. Check coordinates, centre and radius.');
end
if any(~validCoordinates)
    warning('Excluding %d channels with missing scalp coordinates from ROI selection.', ...
        sum(~validCoordinates));
end
roiLabels = string({chanlocs(roiIdx).labels});
fprintf('\nPosterior spatial ROI: %d electrodes, centre [%g, %g], radius %g\n', ...
    numel(roiIdx),roiCentre(1),roiCentre(2),roiRadius);
disp(roiLabels);

% Average selected electrodes within each subject, then average subjects.
roiW = mean(W(:,:,roiIdx),3);
meanROI = mean(roiW,1);
semROI = std(roiW,0,1) ./ sqrt(nIncluded);
n400Window = t >= 300 & t <= 500;
subjectROI300to500 = mean(roiW(:,n400Window),2);

roiFig = figure('Name','Full-trial Word Onset TRF - Posterior Spatial ROI', ...
    'Position',[100,100,1250,500]);
roiLayout = tiledlayout(roiFig,1,2,'TileSpacing','compact');
ax = nexttile(roiLayout);
hold(ax,'on');
lagRow = t(:)';
semPatch = fill(ax,[lagRow,fliplr(lagRow)], ...
    [meanROI+semROI,fliplr(meanROI-semROI)],[0.75,0.85,1], ...
    'EdgeColor','none','FaceAlpha',0.5);
meanLine = plot(ax,t,meanROI,'LineWidth',2);
xline(ax,300,'k--','HandleVisibility','off');
xline(ax,500,'k--','HandleVisibility','off');
yline(ax,0,'k:','HandleVisibility','off');
xlabel(ax,'Latency after word onset (ms)');
ylabel(ax,'Signed TRF weight');
title(ax,sprintf('Posterior spatial ROI (%d electrodes)',numel(roiIdx)));
legend(ax,[meanLine,semPatch],{'Mean TRF','Mean +/- SEM'},'Location','best');
grid(ax,'on');

% Show the selected region separately, without changing the scalp topography.
ax = nexttile(roiLayout);
hold(ax,'on');
circleAngle = linspace(0,2*pi,200);
plot(ax,0.5*cos(circleAngle),0.5*sin(circleAngle),'k-');
plot(ax,[-0.08,0,0.08],[0.49,0.57,0.49],'k-');
scatter(ax,scalpXY(validCoordinates,1),scalpXY(validCoordinates,2),15,[0.6,0.6,0.6],'filled');
scatter(ax,scalpXY(roiIdx,1),scalpXY(roiIdx,2),35,[0,0.447,0.741],'filled');
plot(ax,roiCentre(1)+roiRadius*cos(circleAngle), ...
    roiCentre(2)+roiRadius*sin(circleAngle),'b--');
axis(ax,'equal');
axis(ax,'off');
title(ax,'Selected electrodes');
savefig(roiFig,fullfile(resultDir,'TRF_posteriorSpatialROI.fig'));

%% ============ N400 TOPOGRAPHY ============

% Mean signed TRF weights between 300 and 500 ms
topographyValues = mean(avgTRF(n400Window,:),1);

figure('Name','N400 scalp topography');

topoplot(topographyValues,chanlocs,'electrodes','on');

colorbar;
title('Word Onset TRF - 300 to 500 ms');

savefig(gcf,fullfile(resultDir,'N400_topography.fig'));

%% ============ SAVE GROUP ANALYSIS ============

save(fullfile(resultDir,'N400_groupAnalysis.mat'), ...
    'W','avgTRF','t', ...
    'roiCentre','roiRadius','scalpXY','roiIdx','roiLabels', ...
    'roiW','meanROI','semROI','subjectROI300to500', ...
    'topographyValues','n400Window');

fprintf('\nFull-trial analysis completed.\n');
fprintf('Results saved in:\n%s\n',resultDir);
