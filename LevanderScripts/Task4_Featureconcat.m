%% TASK 4: Feature Concatenation using cndFeatureConcat
% This script demonstrates the use of cndFeatureConcat to combine
% multiple stimulus features into a single concatenated feature.
%
% The function cndFeatureConcat allows flexible grouping of features,
% which is useful when you want to combine related features (e.g., envelope
% and envelope derivative as acoustic features).

clear
close all
clc

%% 1. Setup - Add required paths
disp('=== TASK 4: Feature Concatenation using cndFeatureConcat ===');
disp(' ');

addpath('..\libs\cnsp_utils');
addpath('..\libs\cnsp_utils\cnd');
addpath('..\libs\mTRF-Toolbox_v2\mtrf');
addpath('..\libs\NoiseTools');
addpath('..\libs\eeglab');

% Initialize EEGLAB
try
    eeglab nogui;
    close
catch
    disp('EEGLAB initialization warning - continuing anyway');
end

%% 2. Load data for examples
disp('Loading stimulus data...');

% Try to load from current directory
if exist('dataStim.mat', 'file')
    load('dataStim.mat', 'stim');
    disp('✓ Loaded dataStim.mat from current directory');
else
    % Try alternative path
    dataMainFolder = '../datasets/LalorNatSpeech/dataCND/';
    stimFilename = fullfile(dataMainFolder, 'dataStim.mat');
    if exist(stimFilename, 'file')
        load(stimFilename, 'stim');
        disp(['✓ Loaded: ', stimFilename]);
    else
        error('Could not find dataStim.mat. Please ensure the file exists.');
    end
end

% Display stimulus structure info
disp(' ');
disp('=== STIMULUS DATA INFO ===');
disp(['Sampling frequency: ', num2str(stim.fs), ' Hz']);
disp(['Number of features: ', num2str(size(stim.data, 1))]);
disp(['Number of trials: ', num2str(size(stim.data, 2))]);

if isfield(stim, 'names')
    disp('Feature names:');
    for i = 1:length(stim.names)
        disp(['  Feature ', num2str(i), ': ', stim.names{i}]);
    end
end

% Display feature dimensions
disp(' ');
disp('Feature dimensions (from first trial):');
for i = 1:size(stim.data, 1)
    nDims = size(stim.data{i, 1}, 2);
    disp(['  Feature ', num2str(i), ': ', num2str(nDims), ' dimension(s)']);
end
disp(' ');

%% 3. VISUALIZATION: Show original data structure
disp('=== VISUALIZATION: Original Data Structure ===');
disp(' ');

% Take first 5 samples from first trial for clear visualization
nShow = 5;
fprintf('Original Data Structure (first %d samples, trial 1):\n', nShow);
fprintf('%s\n', repmat('=', 1, 80));

% Create a table-like display
fprintf('\n%-12s', 'Sample #');
for i = 1:size(stim.data, 1)
    fprintf('| %-20s ', ['Feature ' num2str(i)]);
end
fprintf('\n%s\n', repmat('-', 1, 80));

for s = 1:nShow
    fprintf('%-12d', s);
    for i = 1:size(stim.data, 1)
        data = stim.data{i, 1}(s, :);
        if length(data) == 1
            fprintf('| %8.4f           ', data);
        else
            fprintf('| [');
            for d = 1:length(data)
                fprintf('%8.4f', data(d));
                if d < length(data)
                    fprintf(', ');
                end
            end
            fprintf(']    ');
        end
    end
    fprintf('\n');
end
fprintf('\n');

% Show dimensionality info
fprintf('\nDimensionality Summary:\n');
fprintf('%s\n', repmat('-', 1, 50));
totalCols = 0;
for i = 1:size(stim.data, 1)
    nCols = size(stim.data{i, 1}, 2);
    totalCols = totalCols + nCols;
    fprintf('Feature %d: %s → %d column(s)\n', i, stim.names{i}, nCols);
end
fprintf('Total columns across all features: %d\n', totalCols);
fprintf('\n');

%% 4. EXAMPLE 1: Simple concatenation {[1],[2],[3]}
disp('=== EXAMPLE 1: Simple Concatenation {[1],[2],[3]} ===');
disp(' ');

if size(stim.data, 1) >= 3
    stimIdxs1 = {[1], [2], [3]};
    [stimFeature1, feaGroupIdxs1] = cndFeatureConcat(stim, stimIdxs1);
    
    fprintf('Concatenation groups:\n');
    fprintf('  Group 1: Feature 1\n');
    fprintf('  Group 2: Feature 2\n');
    fprintf('  Group 3: Feature 3\n\n');
    
    fprintf('Resulting data (first %d samples):\n', nShow);
    fprintf('%s\n', repmat('=', 1, 80));
    
    % Create header with group boundaries
    fprintf('%-10s', 'Sample #');
    cumCols = 0;
    for g = 1:3
        nCols = size(stim.data{g, 1}, 2);
        for c = 1:nCols
            fprintf('|  G%d-D%d  ', g, c);
        end
        cumCols = cumCols + nCols;
        if g < 3
            fprintf('|');
        end
    end
    fprintf('\n%s\n', repmat('-', 1, 80));
    
    % Print data
    for s = 1:nShow
        fprintf('%-10d', s);
        col = 1;
        for g = 1:3
            nCols = size(stim.data{g, 1}, 2);
            for c = 1:nCols
                fprintf('|  %8.4f', stimFeature1.data{1}(s, col));
                col = col + 1;
            end
            if g < 3
                fprintf('|');
            end
        end
        fprintf('\n');
    end
    
    fprintf('\nGroup indices per dimension:\n');
    fprintf('  [');
    for i = 1:min(20, length(feaGroupIdxs1))
        if i > 1, fprintf(' '); end
        fprintf('%d', feaGroupIdxs1(i));
        if mod(i, 5) == 0 && i < length(feaGroupIdxs1)
            fprintf('\n   ');
        end
    end
    if length(feaGroupIdxs1) > 20
        fprintf(' ...');
    end
    fprintf(']\n');
    fprintf('  (Each number shows which group each dimension belongs to)\n\n');
    
    % Show concatenation visually
    figure('Position', [100, 100, 800, 400]);
    
    subplot(1,2,1);
    % Show original as blocks
    imagesc([stim.data{1,1}(1:20,1) * 2, stim.data{2,1}(1:20,1) * 2 + 1, stim.data{3,1}(1:20,1) * 2 + 2]);
    title('Original Features (color-coded)');
    xlabel('Feature Dimension');
    ylabel('Sample');
    colormap(parula);
    colorbar;
    
    subplot(1,2,2);
    imagesc(stimFeature1.data{1}(1:20, :));
    title('Concatenated Feature');
    xlabel('Feature Dimension');
    ylabel('Sample');
    colorbar;
    
    sgtitle('Feature Concatenation: Simple Grouping {[1],[2],[3]}');
end

%% 5. EXAMPLE 2: Grouped concatenation {[1,2],[3]}
disp(' ');
disp('=== EXAMPLE 2: Grouped Concatenation {[1,2],[3]} ===');
disp(' ');

if size(stim.data, 1) >= 3
    stimIdxs2 = {[1, 2], [3]};
    [stimFeature2, feaGroupIdxs2] = cndFeatureConcat(stim, stimIdxs2);
    
    fprintf('Concatenation groups:\n');
    fprintf('  Group 1: Features 1 and 2 together\n');
    fprintf('  Group 2: Feature 3\n\n');
    
    fprintf('Resulting data (first %d samples):\n', nShow);
    fprintf('%s\n', repmat('=', 1, 80));
    
    % Create header with group boundaries
    fprintf('%-10s', 'Sample #');
    cumCols = 0;
    
    % Group 1 (Features 1 & 2)
    nCols1 = size(stim.data{1,1},2) + size(stim.data{2,1},2);
    for c = 1:nCols1
        fprintf('|  G1-D%-2d ', c);
    end
    fprintf('|');
    
    % Group 2 (Feature 3)
    nCols2 = size(stim.data{3,1},2);
    for c = 1:nCols2
        fprintf('  G2-D%-2d |', c);
    end
    fprintf('\n%s\n', repmat('-', 1, 100));
    
    % Print data
    for s = 1:nShow
        fprintf('%-10d', s);
        col = 1;
        
        % Group 1 data
        for c = 1:nCols1
            fprintf('|   %8.4f', stimFeature2.data{1}(s, col));
            col = col + 1;
        end
        fprintf('|');
        
        % Group 2 data
        for c = 1:nCols2
            fprintf('   %8.4f |', stimFeature2.data{1}(s, col));
            col = col + 1;
        end
        fprintf('\n');
    end
    
    fprintf('\nGroup indices per dimension:\n');
    fprintf('  [');
    for i = 1:min(20, length(feaGroupIdxs2))
        if i > 1, fprintf(' '); end
        fprintf('%d', feaGroupIdxs2(i));
        if mod(i, 5) == 0 && i < length(feaGroupIdxs2)
            fprintf('\n   ');
        end
    end
    if length(feaGroupIdxs2) > 20
        fprintf(' ...');
    end
    fprintf(']\n');
    fprintf('  (Notice: First %d dimensions belong to Group 1, rest to Group 2)\n\n', nCols1);
    
    % Show that data is identical to Example 1
    if exist('stimFeature1', 'var')
        dataMatch = isequal(stimFeature1.data{1}, stimFeature2.data{1});
        if dataMatch
            fprintf('✓ VERIFICATION: Data is IDENTICAL to Example 1\n');
            fprintf('  (Only the grouping information differs)\n\n');
        end
    end
    
    % Visual comparison
    figure('Position', [100, 100, 1000, 300]);
    
    subplot(1,3,1);
    bar(feaGroupIdxs1(1:min(30, length(feaGroupIdxs1))), 'FaceColor', [0.2 0.4 0.8]);
    title('Group Indices - Example 1');
    xlabel('Dimension');
    ylabel('Group');
    ylim([0 4]);
    
    subplot(1,3,2);
    bar(feaGroupIdxs2(1:min(30, length(feaGroupIdxs2))), 'FaceColor', [0.8 0.4 0.2]);
    title('Group Indices - Example 2');
    xlabel('Dimension');
    ylabel('Group');
    ylim([0 4]);
    
    subplot(1,3,3);
    % Show first few samples as a heatmap
    imagesc(stimFeature2.data{1}(1:10, 1:min(20, size(stimFeature2.data{1},2))));
    title('Concatenated Data (first 10 samples)');
    xlabel('Dimension');
    ylabel('Sample');
    colorbar;
    
    sgtitle('Grouping Comparison: Same Data, Different Group Indices');
end

%% 6. SUMMARY VISUALIZATION
disp(' ');
disp('=== SUMMARY: How Concatenation Works ===');
disp(' ');

% Create a simple schematic figure
figure('Position', [100, 100, 1200, 400]);

% Show original features as separate matrices
subplot(1,3,1);
hold on;
rectangle('Position', [1, 1, 2, 5], 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'b', 'LineWidth', 2);
text(1.5, 5.5, 'Feature 1', 'HorizontalAlignment', 'center');
text(1.5, 3, '2 cols', 'HorizontalAlignment', 'center');
rectangle('Position', [4, 1, 1, 5], 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'r', 'LineWidth', 2);
text(4.5, 5.5, 'Feature 2', 'HorizontalAlignment', 'center');
text(4.5, 3, '1 col', 'HorizontalAlignment', 'center');
rectangle('Position', [6, 1, 3, 5], 'FaceColor', [0.8 1 0.8], 'EdgeColor', 'g', 'LineWidth', 2);
text(7.5, 5.5, 'Feature 3', 'HorizontalAlignment', 'center');
text(7.5, 3, '3 cols', 'HorizontalAlignment', 'center');
xlim([0 10]);
ylim([0 7]);
title('Original Features');
axis off;

% Show concatenation process
subplot(1,3,2);
hold on;
% Arrows
annotation('arrow', [0.4 0.45], [0.6 0.6], 'LineWidth', 2);
annotation('arrow', [0.4 0.45], [0.5 0.5], 'LineWidth', 2);
annotation('arrow', [0.4 0.45], [0.4 0.4], 'LineWidth', 2);

rectangle('Position', [1, 2, 2, 2], 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'b');
text(2, 4.2, 'F1', 'HorizontalAlignment', 'center');
rectangle('Position', [1, 4, 2, 2], 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'b');
rectangle('Position', [4, 2, 1, 4], 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'r');
text(4.5, 4.5, 'F2', 'HorizontalAlignment', 'center', 'Rotation', 90);
rectangle('Position', [6, 1, 3, 6], 'FaceColor', [0.8 1 0.8], 'EdgeColor', 'g');
text(7.5, 7.2, 'F3', 'HorizontalAlignment', 'center');
xlim([0 10]);
ylim([0 8]);
title('Stacking Features');
axis off;

% Show concatenated result
subplot(1,3,3);
% Draw concatenated matrix
hold on;
% Feature 1 part
rectangle('Position', [1, 1, 2, 5], 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'b', 'LineWidth', 2);
text(2, 5.5, 'F1', 'HorizontalAlignment', 'center');
% Feature 2 part
rectangle('Position', [3, 1, 1, 5], 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'r', 'LineWidth', 2);
text(3.5, 5.5, 'F2', 'HorizontalAlignment', 'center');
% Feature 3 part
rectangle('Position', [4, 1, 3, 5], 'FaceColor', [0.8 1 0.8], 'EdgeColor', 'g', 'LineWidth', 2);
text(5.5, 5.5, 'F3', 'HorizontalAlignment', 'center');

% Add dimension markers
for i = 1:6
    text(i+0.5, 0.5, sprintf('d%d', i), 'HorizontalAlignment', 'center', 'FontSize', 8);
end
xlim([0 8]);
ylim([0 7]);
title('Concatenated Feature');
axis off;

sgtitle('Feature Concatenation: Side-by-Side Stacking');

%% 7. Key Takeaways
disp(' ');
disp('=== KEY TAKEAWAYS ===');
disp(' ');
disp('  • cndFeatureConcat stacks feature columns side-by-side');
disp('  • Different groupings produce IDENTICAL data arrays');
disp('  • The group indices output tracks which original group each column came from');
disp('  • This is useful for interpreting model weights later');
disp('  • All features must have the same number of time samples');
disp(' ');

disp('=== TASK 4 COMPLETE ===');

%% QUICK DIAGNOSTIC: Show non-zero values before and after concatenation
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('DATA VERIFICATION: Showing non-zero values in scientific notation\n');
fprintf('%s\n', repmat('=', 1, 80));

%% BEFORE CONCATENATION - Show actual values (not rounded)
fprintf('\n▶ ORIGINAL DATA (samples 1-5, scientific notation):\n');
fprintf('%s\n', repmat('-', 1, 60));

for sample = 1:5
    fprintf('Sample %d:\n', sample);
    for feat = 1:4
        data = stim.data{feat, 1}(sample, :);
        if feat == 1
            fprintf('  Feature 1 (envelope): %.6e\n', data(1));
        elseif feat == 2
            fprintf('  Feature 2 (onsets):   %.6e\n', data(1));
        elseif feat == 3
            fprintf('  Feature 3 (spec):     ');
            for d = 1:min(5, length(data))
                fprintf('%.2e  ', data(d));
            end
            if length(data) > 5, fprintf('...'); end
            fprintf('\n');
        elseif feat == 4
            fprintf('  Feature 4 (phonetic): ');
            for d = 1:min(5, length(data))
                fprintf('%.2e  ', data(d));
            end
            if length(data) > 5, fprintf('...'); end
            fprintf('\n');
        end
    end
    fprintf('\n');
end

%% CONCATENATE
fprintf('\n▶ AFTER CONCATENATION (samples 1-5, scientific notation):\n');
fprintf('%s\n', repmat('-', 1, 80));

% Example 1: Simple concatenation
stimIdxs1 = {[1], [2], [3]};
[stimFeature1, feaGroupIdxs1] = cndFeatureConcat(stim, stimIdxs1);

% Show dimensions
fprintf('\nConcatenated matrix dimensions: [%d samples × %d columns]\n', ...
    size(stimFeature1.data{1}, 1), size(stimFeature1.data{1}, 2));
fprintf('Group mapping: ');
for g = 1:3
    nCols = sum(feaGroupIdxs1 == g);
    fprintf('Group%d:%d cols  ', g, nCols);
end
fprintf('\n\n');

% Show first 5 samples with grouping
concatData = stimFeature1.data{1}(1:5, :);
fprintf('Sample  | ');
for col = 1:size(concatData, 2)
    fprintf('Col%d(G%d)  ', col, feaGroupIdxs1(col));
    if mod(col, 5) == 0 && col < size(concatData, 2)
        fprintf('\n        | ');
    end
end
fprintf('\n%s\n', repmat('-', 1, 80));

for sample = 1:5
    fprintf('  %2d    | ', sample);
    for col = 1:size(concatData, 2)
        fprintf('%.2e  ', concatData(sample, col));
        if mod(col, 5) == 0 && col < size(concatData, 2)
            fprintf('\n        | ');
        end
    end
    fprintf('\n');
end

%% VERIFICATION: Show summary statistics
fprintf('\n▶ SUMMARY STATISTICS (full trial):\n');
fprintf('%s\n', repmat('-', 1, 60));

fprintf('Before concatenation:\n');
for feat = 1:4
    data = stim.data{feat, 1}(:);
    fprintf('  Feature %d: range=[%.2e, %.2e], mean=%.2e, std=%.2e\n', ...
        feat, min(data), max(data), mean(data), std(data));
end

fprintf('\nAfter concatenation:\n');
data = stimFeature1.data{1}(:);
fprintf('  Combined: range=[%.2e, %.2e], mean=%.2e, std=%.2e\n', ...
    min(data), max(data), mean(data), std(data));

%% SHOW that Example 1 and 2 are identical
fprintf('\n▶ VERIFICATION: Example 1 vs Example 2\n');
fprintf('%s\n', repmat('-', 1, 60));

stimIdxs2 = {[1,2], [3]};
[stimFeature2, feaGroupIdxs2] = cndFeatureConcat(stim, stimIdxs2);

data1 = stimFeature1.data{1};
data2 = stimFeature2.data{1};

if isequal(data1, data2)
    fprintf('✓ Example 1 and Example 2 produce IDENTICAL data\n');
    fprintf('  Max difference: %.2e\n', max(abs(data1(:) - data2(:))));
else
    fprintf('✗ Data differs! Check implementation\n');
end

%% VISUALIZE the non-zero pattern
figure('Position', [100, 100, 1000, 400]);

subplot(1,3,1);
imagesc(stim.data{1,1}(1:100,1)');
title('Feature 1: Envelope');
xlabel('Sample'); ylabel('Amplitude');
colorbar;

subplot(1,3,2);
imagesc(stim.data{3,1}(1:100,1:8)');
title('Feature 3: Spectrogram (first 8 bands)');
xlabel('Sample'); ylabel('Frequency band');
colorbar;

subplot(1,3,3);
imagesc(stimFeature1.data{1}(1:100, 1:10)');
title('Concatenated (first 10 cols)');
xlabel('Sample'); ylabel('Column');
colorbar;

sgtitle('Non-zero values visualization');