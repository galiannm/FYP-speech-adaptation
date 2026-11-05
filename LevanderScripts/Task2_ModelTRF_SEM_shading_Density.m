%% TASK 2: TRF Weights Plotting
clear; close all; clc;

%% 1. Setup
disp('=== TASK 2: TRF Weights Plotting ===');
addpath('..\libs\cnsp_utils');
addpath('..\libs\cnsp_utils\cnd');
addpath('..\libs\mTRF-Toolbox_v2\mtrf');
addpath('..\libs\NoiseTools');
addpath('..\libs\eeglab');

try
    eeglab nogui; close;
catch
    disp('EEGLAB initialization warning');
end

%% 2. Load TRF model data
disp('Loading TRF model data...');
matFile = fullfile('resultFigures', '1-8Hz', 'TRF_nSub19_condIdx1_stimIdx1.mat');
if ~exist(matFile, 'file')
    matFile = "C:\Users\Dariu\Downloads\CNSP-resources-main\CNSP-resources-main\CNSP\standardAnalysisScripts\resultFigures\1-8Hz\TRF_nSub19_condIdx1_stimIdx1.mat";
end

if exist(matFile, 'file')
    load(matFile, 'modelAll', 'chanlocs');
    disp(['✓ Loaded: ', num2str(length(modelAll)), ' subjects']);
else
    warning('Creating dummy data');
    [modelAll, chanlocs] = createDummyTRFData();
end

%% 3. Get channel names
allLabels = {chanlocs.labels};
disp('Available channels:');
disp(allLabels(1:min(10, length(allLabels))));

% Select example channels
selectedChannels = {};
if any(ismember(allLabels, {'A1', 'A2', 'C1', 'C2'}))
    selectedChannels = {'A1', 'A2'};
    disp('Using A1, A2 (like in screenshot)');
elseif any(ismember(allLabels, {'Cz', 'Pz', 'Fz'}))
    selectedChannels = {'Cz', 'Pz'};
    disp('Using Cz, Pz');
elseif length(allLabels) >= 2
    selectedChannels = allLabels(1:2);
    disp(['Using first two channels: ', strjoin(selectedChannels, ', ')]);
end

disp(' ');

%% 4. DEMONSTRATION - ONLY 5 CLEAR FIGURES

% Figure 1: Two channels with SEM
disp('Creating Figure 1: Two channels with SEM shading');
figure('Color', 'w', 'Position', [100 100 800 600], 'Name', 'Two Channels with SEM');
plotTRFweightsSEM(modelAll, chanlocs, 1, selectedChannels, 1, 18);

% Figure 2: Single channel with SEM
disp('Creating Figure 2: Single channel with SEM');
figure('Color', 'w', 'Position', [100 100 600 500], 'Name', 'Single Channel');
plotTRFweightsSEM(modelAll, chanlocs, 1, selectedChannels(1), 1, 18);

% Figure 3: Two channels with density lines
disp('Creating Figure 3: Two channels with density lines');
figure('Color', 'w', 'Position', [100 100 800 600], 'Name', 'Density Lines');
plotTRFweightsSEM(modelAll, chanlocs, 2, selectedChannels, 1, 18, 3, 'spread');

% Figure 4: All channels simple mean
disp('Creating Figure 4: All channels (mean only)');
figure('Color', 'w', 'Position', [100 100 800 600], 'Name', 'All Channels');
plotTRFweightsSEM(modelAll, chanlocs, 0, [], 1, 18);

% Figure 5: Comparison SEM vs Density
disp('Creating Figure 5: SEM vs Density comparison');
figure('Color', 'w', 'Position', [100 100 1200 500], 'Name', 'Comparison');

% Left subplot: SEM
subplot(1, 2, 1);
plotTRFweightsSEM(modelAll, chanlocs, 1, selectedChannels, 1, 16);

% Right subplot: Density
subplot(1, 2, 2);
plotTRFweightsSEM(modelAll, chanlocs, 2, selectedChannels, 1, 16, 3, 'spread');

% Add overall title
sgtitle('Visualization Comparison', 'FontSize', 16, 'FontWeight', 'bold');

disp(' ');
disp('=== TASK 2 COMPLETE - 5 Figures Generated ===');

%% DUMMY DATA GENERATOR (keep at bottom)
function [modelAll, chanlocs] = createDummyTRFData()
    nSubjects = 19;
    nChannels = 32;
    nTimePoints = 200;
    time = linspace(-100, 400, nTimePoints);
    
    % Create channel labels
    chanlocs = struct();
    letters = {'A', 'B', 'C', 'D'};
    idx = 1;
    for l = 1:length(letters)
        for n = 1:8
            chanlocs(idx).labels = sprintf('%s%d', letters{l}, n);
            idx = idx + 1;
            if idx > nChannels, break; end
        end
        if idx > nChannels, break; end
    end
    
    % Create TRF data
    modelAll = struct();
    for s = 1:nSubjects
        weights = zeros(nChannels, nTimePoints);
        for ch = 1:nChannels
            if contains(chanlocs(ch).labels, 'A')
                peak1 = 80 + randn*8;
                peak2 = 180 + randn*15;
                amp1 = 1.0 + randn*0.2;
                amp2 = -0.6 + randn*0.15;
            elseif contains(chanlocs(ch).labels, 'C')
                peak1 = 120 + randn*10;
                peak2 = 220 + randn*20;
                amp1 = 0.8 + randn*0.15;
                amp2 = -0.4 + randn*0.1;
            else
                peak1 = 100 + randn*10;
                peak2 = 200 + randn*15;
                amp1 = 0.9 + randn*0.18;
                amp2 = -0.5 + randn*0.12;
            end
            
            waveform = amp1 * exp(-(time-peak1).^2/(2*25^2)) + ...
                       amp2 * exp(-(time-peak2).^2/(2*30^2));
            
            noiseLevel = 0.07 + randn*0.02;
            weights(ch,:) = waveform + randn(1, nTimePoints)*noiseLevel;
        end
        modelAll(s).w = weights;
        modelAll(s).t = time;
    end
    disp('Created dummy TRF data');
end