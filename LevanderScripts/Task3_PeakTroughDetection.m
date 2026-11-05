%% TRF Peaks Analysis - Example Usage Script
clear; close all; clc;

%% Add paths and load data
% run addPathDependencies.m  % Uncomment if needed

matFile = fullfile('resultFigures','1-8Hz','TRF_nSub19_condIdx1_stimIdx1.mat');
fprintf('Loading: %s\n', matFile);
load(matFile, 'modelAll', 'chanlocs');
fprintf('Loaded %d models with %d channels\n', length(modelAll), length(chanlocs));

%% ========================================================================
%% EXAMPLE 1: Basic usage - Get peaks and interactive plot
%% ========================================================================
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('EXAMPLE 1: Basic interactive plot\n');
fprintf('%s\n', repmat('=', 1, 70));

[peaks1, fig1] = cnsp_plotTRFPeaks(modelAll, chanlocs, ...
    'nPeaks', 3, ...
    'nTroughs', 2, ...
    'timeRange', [-100 600], ...
    'figTitle', 'TRF Analysis: 3 Peaks, 2 Troughs');

%% ========================================================================
%% EXAMPLE 2: Data-only mode - Just get peak data without plot
%% ========================================================================
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('EXAMPLE 2: Data-only mode (no plot)\n');
fprintf('%s\n', repmat('=', 1, 70));

peaks_data = cnsp_plotTRFPeaks(modelAll, chanlocs, ...
    'displayTopoplot', false, ...  % Don't show plot
    'nPeaks', 3, ...
    'nTroughs', 2, ...
    'timeRange', [-100 600]);

% Display results
fprintf('\n--- Detected Peaks ---\n');
for i = 1:peaks_data.nPeaks
    fprintf('  P%d: %.1f ms, value=%.3f\n', ...
        i, peaks_data.positive.times(i), peaks_data.positive.values(i));
end

fprintf('\n--- Detected Troughs ---\n');
for i = 1:peaks_data.nTroughs
    fprintf('  N%d: %.1f ms, value=%.3f\n', ...
        i, peaks_data.negative.times(i), peaks_data.negative.values(i));
end

%% ========================================================================
%% EXAMPLE 3: Detect all peaks and troughs
%% ========================================================================
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('EXAMPLE 3: Detect all peaks and troughs\n');
fprintf('%s\n', repmat('=', 1, 70));

[peaks_all, fig2] = cnsp_plotTRFPeaks(modelAll, chanlocs, ...
    'nPeaks', 10, ...      % Large number to get all peaks
    'nTroughs', 10, ...    % Large number to get all troughs
    'figTitle', 'TRF Analysis: All Peaks/Troughs');

fprintf('\nFound %d peaks and %d troughs\n', ...
    peaks_all.nPeaks, peaks_all.nTroughs);

%% ========================================================================
%% EXAMPLE 4: Access topography data for specific time points
%% ========================================================================
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('EXAMPLE 4: Access topography data\n');
fprintf('%s\n', repmat('=', 1, 70));

% Get data only
peaks = cnsp_plotTRFPeaks(modelAll, chanlocs, ...
    'displayTopoplot', false, ...
    'nPeaks', 3, ...
    'nTroughs', 2);

% Access individual topographies
for i = 1:peaks.nTotal
    fprintf('\nTopography for %s (%.1f ms):\n', ...
        peaks.labels{i}, peaks.timepoints(i));
    topo = peaks.topographies(i, :);
    fprintf('  Min: %.3f, Max: %.3f, Mean: %.3f\n', ...
        min(topo), max(topo), mean(topo));
end

%% ========================================================================
%% Summary
%% ========================================================================
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 70));

fprintf('\nInteractive features:\n');
fprintf('  • Click ANYWHERE on TRF line → shows topography at that time\n');
fprintf('  • Hover mouse → shows time indicator line\n');
fprintf('  • Detected peaks (P#) and troughs (N#) are marked\n');
fprintf('  • Click on peaks/troughs to see their topographies\n');

fprintf('\n=== Done ===\n');