%% TASK 6: Quick Test of cndCheck Function
% Simple test script for the cndCheck function

clear; close all; clc;

%% 1. Setup - Add required paths
disp('TASK 6: Testing cndCheck Function');
disp('=================================');
fprintf('\n');

addpath('..\libs\cnsp_utils');
addpath('..\libs\cnsp_utils\cnd');

%% 2. Define test files
dataFolder = '../datasets/LalorNatSpeech/dataCND/';
stimFile = fullfile(dataFolder, 'dataStim.mat');
neuralFile = fullfile(dataFolder, 'dataSub1.mat');

%% 3. Test 1: Basic validation with verbose output
disp('TEST 1: Basic validation (verbose mode)');
disp('---------------------------------------');

if exist(stimFile, 'file') && exist(neuralFile, 'file')
    [summary1, results1, passed1] = cndCheck(stimFile, neuralFile, 1);
    
    fprintf('\nResults:\n');
    fprintf('  Overall passed: %s\n', string(passed1));
    fprintf('  Stimulus valid: %s\n', string(results1.stimValid));
    fprintf('  Neural valid: %s\n', string(results1.neuralValid));
else
    fprintf('❌ Files not found\n');
end

%% 4. Test 2: Quick check without verbose output
disp('\n\nTEST 2: Quick validation (silent mode)');
disp('---------------------------------------');

if exist(stimFile, 'file') && exist(neuralFile, 'file')
    [summary2, results2, passed2] = cndCheck(stimFile, neuralFile, 0);
    
    fprintf('Summary table:\n');
    disp(summary2);
    fprintf('Overall status: %s\n', string(passed2));
else
    fprintf('❌ Files not found\n');
end

%% 5. Test 3: Validate only stimulus data
disp('\n\nTEST 3: Stimulus-only validation');
disp('----------------------------------');

if exist(stimFile, 'file')
    [summary3, results3, passed3] = cndCheck(stimFile, [], 1);
    fprintf('Stimulus validation: %s\n', string(results3.stimValid));
else
    fprintf('❌ Stimulus file not found\n');
end

%% 6. Test 4: Validate only neural data
disp('\n\nTEST 4: Neural-only validation');
disp('--------------------------------');

if exist(neuralFile, 'file')
    [summary4, results4, passed4] = cndCheck([], neuralFile, 1);
    fprintf('Neural validation: %s\n', string(results4.neuralValid));
else
    fprintf('❌ Neural file not found\n');
end

%% 7. Test 5: Test with custom options
disp('\n\nTEST 5: Custom options test');
disp('----------------------------');

if exist(stimFile, 'file') && exist(neuralFile, 'file')
    options = struct();
    options.tolerance = 0.02; % 2% tolerance for trial lengths
    options.maxSamples = 3;   % Check only 3 samples
    
    [summary5, results5, passed5] = cndCheck(stimFile, neuralFile, 1, options);
    fprintf('Custom options test completed: %s\n', string(passed5));
else
    fprintf('❌ Files not found\n');
end

%% 8. Summary
fprintf('\n\n=== TEST SUMMARY ===\n');
fprintf('All tests completed.\n');
fprintf('Check above for any errors or warnings.\n');