clear;
clc;

%% Load stimulus
load('dataStim.mat');

%% Load one participant
load('dataSub1.mat');

%% Inspect structures
disp(stim)
disp(eeg)

%% Basic information
fprintf('Stimulus sampling rate: %g Hz\n', stim.fs);
fprintf('EEG sampling rate: %g Hz\n', eeg.fs);
fprintf('Number of EEG trials: %d\n', length(eeg.data));
fprintf('Number of stimulus trials: %d\n', size(stim.data,2));

%% Look at dimensions of trial 1
size(stim.data{1,1})
size(stim.data{2,1})
size(stim.data{3,1})
size(eeg.data{1})