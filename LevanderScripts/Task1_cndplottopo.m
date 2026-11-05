clear
close all
clc

%% 1. Setup and load data
disp('=== EEG Topography Plotting Examples ===');
disp(' ');
run addPathDependencies.m
% Add necessary paths
addpath('..\libs\cnsp_utils');
addpath('..\libs\cnsp_utils\cnd');
addpath('..\libs\mTRF-Toolbox_v2\mtrf');
addpath('..\libs\NoiseTools');
addpath('..\libs\eeglab');


datasetName = 'LalorNatSpeech';
freqBandName = '1-8Hz'; 
resultsPath = ['./resultFigures/', freqBandName, '/'];

% 1. Load the pre-calculated results
resFile = dir([resultsPath, 'TRF_nSub*_condIdx1_stimIdx1.mat']);
if isempty(resFile)
    error('No results file found.');
else
    load(fullfile(resFile.folder, resFile.name), 'rAllElec', 'chanlocs');
end

meanR = mean(rAllElec, 2)'; 

% ---------------------------------------------------------
% OPTION A: Static Mean Plot
% ---------------------------------------------------------
cndplottopo(meanR, chanlocs, 'title', 'Mean R-Values', 'legend', 'R', 'hover', 'on');
caxis([-0.05, 0.05]); 

% ---------------------------------------------------------
% OPTION B: Interactive Subject Browser
% Uses 'frameName' to show "Subject X/N" in the title.
% ---------------------------------------------------------
cndplottopo(rAllElec', chanlocs, ...
    'title', 'Subject Correlation', ...
    'frameName', 'Subject', ...
    'labels', 'on');
caxis([-0.1, 0.1]);

% ---------------------------------------------------------
% OPTION C: Channel Identification Mode
% ---------------------------------------------------------
cndplottopo(meanR, chanlocs, 'title', 'Click to Identify', 'labels', 'both');

% ---------------------------------------------------------
% OPTION D: Blank Topoplot (Labels ON, Hover OFF)
% ---------------------------------------------------------
cndplottopo([], chanlocs, 'title', 'Layout', 'labels', 'on', 'hover', 'off');



