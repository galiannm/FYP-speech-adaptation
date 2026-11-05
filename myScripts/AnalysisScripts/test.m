clear
close all
clc

%% 1. Set paths

projectRoot = '/Users/manon/Documents/Uni/Trinity4/FYP/CNSP-resources/CNSP';

datasetName = 'AAD_KULeuven';

% Add required libraries
addpath(genpath(fullfile(projectRoot, 'libs')));

% Path to preprocessed EEG for participant 1
eegFile = fullfile( ...
    projectRoot, ...
    'datasets2', ...
    datasetName, ...
    'dataCND', ...
    '1-8Hz', ...
    'pre_dataSub1.mat');


%% 2. Load EEG data

S = load(eegFile);
eeg = S.eeg;

fprintf('\n===== EEG INFORMATION =====\n');
fprintf('Dataset: %s\n', datasetName);
fprintf('Sampling rate: %g Hz\n', eeg.fs);
fprintf('Number of trials: %d\n', numel(eeg.data));
fprintf('Number of channels: %d\n', size(eeg.data{1}, 2));


%% 3. Inspect electrode information

fprintf('\n===== CHANNEL LOCATION FIELDS =====\n');
disp(fieldnames(eeg.chanlocs));

% Check that channel labels exist
if ~isfield(eeg.chanlocs, 'labels')
    error('No channel labels found in eeg.chanlocs');
end

% Extract all electrode labels
labels = string({eeg.chanlocs.labels});

fprintf('\n===== ALL ELECTRODE LABELS =====\n');

for i = 1:numel(labels)
    fprintf('%d: %s\n', i, labels(i));
end


%% 4. Identify possible mastoid electrodes

% Common mastoid electrode names
mastoidNames = ["M1", "M2", "A1", "A2", "TP9", "TP10"];

% Find matching labels
matches = ismember(upper(labels), mastoidNames);

fprintf('\n===== POSSIBLE MASTOID CHANNELS =====\n');

if any(matches)
    disp(labels(matches)');
else
    fprintf('No commonly named mastoid channels found.\n');
end


%% 5. Inspect the cndReref function

fprintf('\n===== RE-REFERENCING FUNCTION =====\n');

% Check whether the function is available
funcPath = which('cndReref');

if isempty(funcPath)
    error('cndReref not found. Check your library paths.');
else
    fprintf('Function found at:\n%s\n', funcPath);
end

% Display documentation
help cndReref


%% 6. Test mastoid re-referencing

fprintf('\n===== TESTING MASTOID RE-REFERENCING =====\n');

% Copy EEG data so the original is preserved
eegTest = eeg;

try
    % Apply the same re-referencing used in your TRF script
    eegTest = cndReref(eegTest, 'Mastoids');

    fprintf('Re-referencing completed successfully.\n');

    % Compare channel counts
    fprintf('Channels before: %d\n', size(eeg.data{1}, 2));
    fprintf('Channels after: %d\n', size(eegTest.data{1}, 2));

    % Check whether EEG values changed
    if isequal(eeg.data{1}, eegTest.data{1})
        fprintf('Warning: EEG values are unchanged.\n');
    else
        fprintf('EEG values changed after re-referencing.\n');
    end

catch ME
    fprintf('Re-referencing failed:\n%s\n', ME.message);
end


%% 7. Summary

fprintf('\n===== SUMMARY =====\n');
fprintf('Dataset: %s\n', datasetName);
fprintf('Total EEG channels: %d\n', numel(labels));
fprintf('Possible mastoid channels: %d\n', sum(matches));

if any(matches)
    fprintf('Mastoid candidates: %s\n', ...
        strjoin(labels(matches), ', '));
end

fprintf('\nInspection completed.\n');