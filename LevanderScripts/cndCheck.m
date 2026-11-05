function [summaryTable, checkResults, passedAll] = cndCheck(stimFile, neuralFile, verbose, options)
%CNDCHECK Validate files against CND (Continuous-event Neural Data) format
%   Author: Evander
%
%   [SUMMARYTABLE, CHECKRESULTS, PASSEDALL] = CNDCHECK(STIMFILE, NEURALFILE, VERBOSE, OPTIONS)
%   Validates CND format data files. CND is a standardized format for neural
%   data analysis. Validation covers structure integrity, required fields,
%   data format, internal consistency, and cross-file compatibility.
%
%   INPUTS:
%       stimFile   - Path to stimulus data MAT file or [] (optional)
%       neuralFile - Path to neural data MAT file or [] (optional)
%       verbose    - 1: detailed output, 0: summary only (default=0)
%       options    - Optional configuration structure (see below)
%
%   OUTPUTS:
%       summaryTable - Table with category-wise pass/fail counts
%       checkResults - Structure with detailed validation results
%       passedAll    - Boolean: true if all checks passed
%
%   OPTIONS (optional fields):
%       options.neuralVarName    - Allowed neural variable names
%                                  Default: {'neural','eeg','EEG','erp','neural_data'}
%       options.tolerance        - Trial length tolerance (0-1, default=0.01)
%       options.maxSamples       - Max samples for NaN/Inf checks (default=5)
%       options.requireStimNames - Require stim.names field (default=false)
%       options.requireChanlocs  - Require neural.chanlocs (default=false)
%
%   EXAMPLES:
%       % Basic validation with detailed output
%       [summary, results, passed] = cndCheck('dataStim.mat', 'dataSub1.mat', 1);
%
%       % Quick validation (silent mode)
%       [summary, results, passed] = cndCheck('dataStim.mat', 'dataSub1.mat', 0);
%
%       % Validate stimulus only
%       [summary, results, passed] = cndCheck('dataStim.mat', [], 1);
%
%       % Validate neural only
%       [summary, results, passed] = cndCheck([], 'dataSub1.mat', 1);
%
%       % With custom options
%       options.tolerance = 0.02;
%       options.maxSamples = 3;
%       [summary, results, passed] = cndCheck('stim.mat', 'neural.mat', 1, options);
%
%
%   Version: 1.0.0
%   Last modified: 2024-01-15

% =========================================================================
% INPUT VALIDATION AND SETUP
% =========================================================================

    % Check minimum number of inputs
    if nargin < 1
        error('cndCheck:missingInput', 'At least one input is required');
    end
    
    % Set defaults for missing inputs
    if nargin < 2
        neuralFile = [];
    end
    
    if nargin < 3
        verbose = 0;
    end
    
    if nargin < 4
        options = struct();
    end
    
    % Validate input types
    validateInputs(stimFile, neuralFile, verbose);
    
    % Set default options
    options = setDefaultOptions(options);
    
    % Initialize results structure
    [checkResults, categories] = initializeResults();
    
    % Define all checks
    allChecks = defineChecks(options);
    checkResults.allChecks = allChecks;
    
% =========================================================================
% LOAD DATA FROM FILES
% =========================================================================

    [stim, neural, loadMessages] = loadDataFiles(stimFile, neuralFile, options, verbose);
    checkResults.loadMessages = loadMessages;
    
    % Get data statistics for reporting
    checkResults.dataStats = getDataStatistics(stim, neural);
    
% =========================================================================
% PERFORM VALIDATION CHECKS
% =========================================================================

    if verbose
        fprintf('\n=== CND VALIDATION CHECKS ===\n\n');
    end
    
    % Perform all checks
    [allChecks, checkResults] = performChecks(allChecks, checkResults, ...
        stim, neural, options, verbose);
    
    % Store updated checks
    checkResults.allChecks = allChecks;
    
% =========================================================================
% CREATE SUMMARY AND FINAL RESULTS
% =========================================================================

    % Create summary table
    summaryTable = createSummaryTable(checkResults);
    
    % Calculate final validation flags
    checkResults = calculateValidationFlags(checkResults);
    passedAll = checkResults.passedAll;
    
    % Display results if verbose
    if verbose
        displayResults(summaryTable, checkResults);
    end
    
end % Main function


% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function validateInputs(stimFile, neuralFile, verbose)
% Validate input parameters
    
    % Check stimFile
    if ~isempty(stimFile)
        if ~(ischar(stimFile) || isstring(stimFile))
            error('cndCheck:invalidInput', ...
                'stimFile must be a string, character array, or empty');
        end
        if ~exist(stimFile, 'file')
            error('cndCheck:fileNotFound', ...
                'Stimulus file not found: %s', stimFile);
        end
    end
    
    % Check neuralFile
    if ~isempty(neuralFile)
        if ~(ischar(neuralFile) || isstring(neuralFile))
            error('cndCheck:invalidInput', ...
                'neuralFile must be a string, character array, or empty');
        end
        if ~exist(neuralFile, 'file')
            error('cndCheck:fileNotFound', ...
                'Neural file not found: %s', neuralFile);
        end
    end
    
    % Check verbose flag
    if ~isscalar(verbose) || ~(isnumeric(verbose) || islogical(verbose))
        error('cndCheck:invalidInput', ...
            'verbose must be a scalar numeric or logical value');
    end
end

function options = setDefaultOptions(options)
% Set default values for options structure
    
    defaultOptions = struct(...
        'neuralVarName', {{'neural', 'eeg', 'EEG', 'erp', 'neural_data'}}, ...
        'checkAllTrials', false, ...
        'tolerance', 0.01, ...
        'maxSamples', 5, ...
        'requireStimNames', false, ...
        'requireChanlocs', false);
    
    % Merge user options with defaults
    optionFields = fieldnames(defaultOptions);
    for i = 1:length(optionFields)
        field = optionFields{i};
        if ~isfield(options, field)
            options.(field) = defaultOptions.(field);
        end
    end
end

function [checkResults, categories] = initializeResults()
% Initialize the results structure
    
    checkResults = struct();
    checkResults.checksPassed = 0;
    checkResults.checksTotal = 0;
    checkResults.categories = struct();
    checkResults.loadMessages = {};
    checkResults.dataStats = struct();
    
    % Category definitions
    categories = {
        'Structure Integrity', 'Required Fields', 'Data Format', ...
        'Data Consistency', 'Compatibility'
    };
    
    % Initialize category counters
    for i = 1:length(categories)
        catName = matlab.lang.makeValidName(categories{i});
        checkResults.categories.(catName) = struct(...
            'passed', 0, 'failed', 0, 'total', 0, ...
            'failures', {{}}, 'warnings', {{}});
    end
end

function allChecks = defineChecks(options)
% Define all validation checks
    
    allChecks = {
        % Structure Integrity checks
        {'Structure Integrity', 'Stimulus is a structure', 'stimIsStruct', []};
        {'Structure Integrity', 'Neural is a structure', 'neuralIsStruct', []};
        
        % Required Fields checks
        {'Required Fields', 'Stimulus has data field', 'stimHasData', []};
        {'Required Fields', 'Stimulus has fs field', 'stimHasFs', []};
        {'Required Fields', 'Neural has data field', 'neuralHasData', []};
        {'Required Fields', 'Neural has fs field', 'neuralHasFs', []};
    };
    
    % Add optional field checks based on options
    if options.requireStimNames
        allChecks{end+1} = {'Required Fields', 'Stimulus has names field', 'stimHasNames', []};
    else
        allChecks{end+1} = {'Data Format', 'Stimulus names field optional', 'stimNamesOptional', []};
    end
    
    if options.requireChanlocs
        allChecks{end+1} = {'Required Fields', 'Neural has chanlocs field', 'neuralHasChanlocs', []};
    else
        allChecks{end+1} = {'Data Format', 'Neural chanlocs optional', 'neuralChanlocsOptional', []};
    end
    
    % Data Format checks (continued)
    dataFormatChecks = {
        {'Data Format', 'Stimulus data is cell array', 'stimDataIsCell', []};
        {'Data Format', 'Neural data is cell array', 'neuralDataIsCell', []};
        {'Data Format', 'Stimulus data not empty', 'stimDataNotEmpty', []};
        {'Data Format', 'Neural data not empty', 'neuralDataNotEmpty', []};
        {'Data Format', 'Stimulus data numeric', 'stimDataNumeric', []};
        {'Data Format', 'Neural data numeric', 'neuralDataNumeric', []};
        {'Data Format', 'Stimulus fs positive number', 'stimFsPositive', []};
        {'Data Format', 'Neural fs positive number', 'neuralFsPositive', []};
    };
    
    allChecks = [allChecks; dataFormatChecks];
    
    % Data Consistency checks
    dataConsistencyChecks = {
        {'Data Consistency', 'Stimulus feature dimensions consistent', 'stimDimsConsistent', []};
        {'Data Consistency', 'Neural channels consistent', 'neuralChansConsistent', []};
        {'Data Consistency', 'Stimulus no NaN values', 'stimNoNaN', []};
        {'Data Consistency', 'Neural no NaN values', 'neuralNoNaN', []};
        {'Data Consistency', 'Stimulus no Inf values', 'stimNoInf', []};
        {'Data Consistency', 'Neural no Inf values', 'neuralNoInf', []};
    };
    
    allChecks = [allChecks; dataConsistencyChecks];
    
    % Compatibility checks
    compatibilityChecks = {
        {'Compatibility', 'Sampling frequencies match', 'fsMatch', []};
        {'Compatibility', 'Number of trials match', 'nTrialsMatch', []};
        {'Compatibility', 'Trial lengths compatible', 'trialLengthsCompat', []};
    };
    
    allChecks = [allChecks; compatibilityChecks];
end

function [stim, neural, loadMessages] = loadDataFiles(stimFile, neuralFile, options, verbose)
% Load data from files with error handling
    
    stim = [];
    neural = [];
    loadMessages = {};
    
    % Load stimulus data
    if ~isempty(stimFile)
        try
            if verbose, fprintf('Loading stimulus data: %s\n', stimFile); end
            load(stimFile, 'stim');
            if ~exist('stim', 'var')
                error('cndCheck:missingVariable', ...
                    'File does not contain ''stim'' variable');
            end
            if verbose, fprintf('  ✓ Stimulus data loaded\n'); end
            loadMessages{end+1} = sprintf('Stimulus loaded from: %s', stimFile);
        catch ME
            loadMessages{end+1} = sprintf('Failed to load stimulus: %s', ME.message);
            if verbose
                fprintf('  ✗ Failed to load stimulus: %s\n', ME.message);
            end
        end
    end
    
    % Load neural data
    if ~isempty(neuralFile)
        try
            if verbose, fprintf('Loading neural data: %s\n', neuralFile); end
            load(neuralFile);
            
            % Check for allowed variable names
            varLoaded = false;
            for v = 1:length(options.neuralVarName)
                varName = options.neuralVarName{v};
                if exist(varName, 'var')
                    neural = eval(varName);
                    varLoaded = true;
                    break;
                end
            end
            
            if ~varLoaded
                error('cndCheck:missingVariable', ...
                    'File must contain one of: %s', ...
                    strjoin(options.neuralVarName, ', '));
            end
            
            if verbose, fprintf('  ✓ Neural data loaded\n'); end
            loadMessages{end+1} = sprintf('Neural loaded from: %s', neuralFile);
        catch ME
            loadMessages{end+1} = sprintf('Failed to load neural: %s', ME.message);
            if verbose
                fprintf('  ✗ Failed to load neural: %s\n', ME.message);
            end
        end
    end
    
    if verbose && (~isempty(stimFile) || ~isempty(neuralFile))
        fprintf('\n');
    end
end

function dataStats = getDataStatistics(stim, neural)
% Calculate and return data statistics
    
    dataStats = struct();
    
    % Stimulus statistics
    if ~isempty(stim) && isstruct(stim)
        dataStats.stim = struct();
        if isfield(stim, 'fs')
            dataStats.stim.fs = stim.fs;
        end
        if isfield(stim, 'data') && iscell(stim.data)
            dataStats.stim.nFeatures = size(stim.data, 1);
            dataStats.stim.nTrials = size(stim.data, 2);
            if dataStats.stim.nTrials > 0 && ~isempty(stim.data{1,1})
                dataStats.stim.trialLength = size(stim.data{1,1}, 1);
                dataStats.stim.nFeaturesPerTrial = size(stim.data{1,1}, 2);
            end
        end
    end
    
    % Neural statistics
    if ~isempty(neural) && isstruct(neural)
        dataStats.neural = struct();
        if isfield(neural, 'fs')
            dataStats.neural.fs = neural.fs;
        end
        if isfield(neural, 'data') && iscell(neural.data)
            dataStats.neural.nTrials = length(neural.data);
            if dataStats.neural.nTrials > 0 && ~isempty(neural.data{1})
                dataStats.neural.trialLength = size(neural.data{1}, 1);
                dataStats.neural.nChannels = size(neural.data{1}, 2);
            end
        end
    end
end

function [allChecks, checkResults] = performChecks(allChecks, checkResults, ...
    stim, neural, options, verbose)
% Perform all validation checks
    
    % Get number of checks
    if iscell(allChecks) && ~isempty(allChecks)
        numChecks = length(allChecks);
    else
        numChecks = 0;
    end
    
    % Perform each check
    for i = 1:numChecks
        currentCheck = allChecks{i};
        
        category = currentCheck{1};
        description = currentCheck{2};
        checkName = currentCheck{3};
        
        % Get category struct name
        catStructName = matlab.lang.makeValidName(category);
        
        % Perform the specific check
        [passed, message, severity] = performSingleCheck(...
            checkName, stim, neural, options);
        
        % Store check result
        allChecks{i}{4} = struct(...
            'passed', passed, ...
            'message', message, ...
            'severity', severity, ...
            'suggestion', getSuggestion(checkName, message));
        
        % Update category counters
        catStruct = checkResults.categories.(catStructName);
        catStruct.total = catStruct.total + 1;
        
        if passed
            catStruct.passed = catStruct.passed + 1;
            checkResults.checksPassed = checkResults.checksPassed + 1;
            if verbose
                fprintf('  ✓ %s\n', description);
            end
        else
            catStruct.failed = catStruct.failed + 1;
            if strcmp(severity, 'error')
                catStruct.failures{end+1} = sprintf('%s: %s', description, message);
            else
                catStruct.warnings{end+1} = sprintf('%s: %s', description, message);
            end
            
            if verbose
                if strcmp(severity, 'error')
                    fprintf('  ✗ %s\n', description);
                else
                    fprintf('  ⚠ %s\n', description);
                end
                if ~isempty(message)
                    fprintf('     Reason: %s\n', message);
                end
                suggestion = getSuggestion(checkName, message);
                if ~isempty(suggestion)
                    fprintf('     Suggestion: %s\n', suggestion);
                end
            end
        end
        
        checkResults.categories.(catStructName) = catStruct;
        checkResults.checksTotal = checkResults.checksTotal + 1;
    end
end

function [passed, message, severity] = performSingleCheck(checkName, stim, neural, options)
% Perform a single validation check
    
    passed = false;
    message = '';
    severity = 'error'; % Default to error
    
    switch checkName
        % Structure Integrity checks
        case 'stimIsStruct'
            if ~isempty(stim)
                passed = isstruct(stim);
                if ~passed
                    message = 'stim is not a structure';
                end
            else
                passed = true;
            end
            
        case 'neuralIsStruct'
            if ~isempty(neural)
                passed = isstruct(neural);
                if ~passed
                    message = 'neural is not a structure';
                end
            else
                passed = true;
            end
            
        % Required Fields checks
        case 'stimHasData'
            if ~isempty(stim)
                passed = isfield(stim, 'data');
                if ~passed
                    message = 'Missing stim.data field';
                end
            else
                passed = true;
            end
            
        case 'stimHasFs'
            if ~isempty(stim)
                passed = isfield(stim, 'fs');
                if ~passed
                    message = 'Missing stim.fs field';
                end
            else
                passed = true;
            end
            
        case 'neuralHasData'
            if ~isempty(neural)
                passed = isfield(neural, 'data');
                if ~passed
                    message = 'Missing neural.data field';
                end
            else
                passed = true;
            end
            
        case 'neuralHasFs'
            if ~isempty(neural)
                passed = isfield(neural, 'fs');
                if ~passed
                    message = 'Missing neural.fs field';
                end
            else
                passed = true;
            end
            
        case 'stimHasNames'
            if ~isempty(stim)
                passed = isfield(stim, 'names');
                if ~passed
                    message = 'stim.names field is required but missing';
                end
            else
                passed = true;
            end
            
        case 'neuralHasChanlocs'
            if ~isempty(neural)
                passed = isfield(neural, 'chanlocs');
                if ~passed
                    message = 'neural.chanlocs field is required but missing';
                end
            else
                passed = true;
            end
            
        % Data Format checks
        case 'stimDataIsCell'
            if ~isempty(stim) && isfield(stim, 'data')
                passed = iscell(stim.data);
                if ~passed
                    message = 'stim.data is not a cell array';
                end
            else
                passed = true;
            end
            
        case 'neuralDataIsCell'
            if ~isempty(neural) && isfield(neural, 'data')
                passed = iscell(neural.data);
                if ~passed
                    message = 'neural.data is not a cell array';
                end
            else
                passed = true;
            end
            
        case 'stimDataNotEmpty'
            if ~isempty(stim) && isfield(stim, 'data')
                passed = ~isempty(stim.data);
                if ~passed
                    message = 'stim.data is empty';
                end
            else
                passed = true;
            end
            
        case 'neuralDataNotEmpty'
            if ~isempty(neural) && isfield(neural, 'data')
                passed = ~isempty(neural.data);
                if ~passed
                    message = 'neural.data is empty';
                end
            else
                passed = true;
            end
            
        case 'stimDataNumeric'
            if ~isempty(stim) && isfield(stim, 'data') && ~isempty(stim.data)
                passed = true;
                % Sample multiple cells for checking
                nSamples = min(options.maxSamples, numel(stim.data));
                for idx = 1:nSamples
                    if ~isnumeric(stim.data{idx})
                        passed = false;
                        message = 'stim.data contains non-numeric values';
                        break;
                    end
                end
            else
                passed = true;
            end
            
        case 'neuralDataNumeric'
            if ~isempty(neural) && isfield(neural, 'data') && ~isempty(neural.data)
                passed = true;
                % Sample multiple trials for checking
                nSamples = min(options.maxSamples, length(neural.data));
                for idx = 1:nSamples
                    if ~isnumeric(neural.data{idx})
                        passed = false;
                        message = 'neural.data contains non-numeric values';
                        break;
                    end
                end
            else
                passed = true;
            end
            
        case 'stimFsPositive'
            if ~isempty(stim) && isfield(stim, 'fs')
                passed = isnumeric(stim.fs) && stim.fs > 0;
                if ~passed
                    message = 'stim.fs is not a positive number';
                end
            else
                passed = true;
            end
            
        case 'neuralFsPositive'
            if ~isempty(neural) && isfield(neural, 'fs')
                passed = isnumeric(neural.fs) && neural.fs > 0;
                if ~passed
                    message = 'neural.fs is not a positive number';
                end
            else
                passed = true;
            end
            
        case 'stimNamesOptional'
            if ~isempty(stim)
                if isfield(stim, 'names')
                    passed = iscell(stim.names);
                    if ~passed
                        message = 'stim.names exists but is not a cell array';
                        severity = 'warning';
                    end
                else
                    passed = true; % Optional field, absence is OK
                end
            else
                passed = true;
            end
            
        case 'neuralChanlocsOptional'
            if ~isempty(neural)
                if isfield(neural, 'chanlocs')
                    passed = isstruct(neural.chanlocs);
                    if ~passed
                        message = 'neural.chanlocs exists but is not a structure';
                        severity = 'warning';
                    end
                else
                    passed = true; % Optional field, absence is OK
                end
            else
                passed = true;
            end
            
        % Data Consistency checks
        case 'stimDimsConsistent'
            if ~isempty(stim) && isfield(stim, 'data') && ~isempty(stim.data)
                passed = checkStimDimensions(stim, options);
                if ~passed
                    message = 'Feature dimensions inconsistent across trials';
                end
            else
                passed = true;
            end
            
        case 'neuralChansConsistent'
            if ~isempty(neural) && isfield(neural, 'data') && ~isempty(neural.data)
                passed = checkNeuralChannels(neural, options);
                if ~passed
                    message = 'Inconsistent number of channels across trials';
                end
            else
                passed = true;
            end
            
        case 'stimNoNaN'
            if ~isempty(stim) && isfield(stim, 'data') && ~isempty(stim.data)
                passed = checkForNaN(stim.data, 'stimulus', options);
                if ~passed
                    message = 'stim.data contains NaN values';
                end
            else
                passed = true;
            end
            
        case 'neuralNoNaN'
            if ~isempty(neural) && isfield(neural, 'data') && ~isempty(neural.data)
                passed = checkForNaN(neural.data, 'neural', options);
                if ~passed
                    message = 'neural.data contains NaN values';
                end
            else
                passed = true;
            end
            
        case 'stimNoInf'
            if ~isempty(stim) && isfield(stim, 'data') && ~isempty(stim.data)
                passed = checkForInf(stim.data, 'stimulus', options);
                if ~passed
                    message = 'stim.data contains Inf values';
                end
            else
                passed = true;
            end
            
        case 'neuralNoInf'
            if ~isempty(neural) && isfield(neural, 'data') && ~isempty(neural.data)
                passed = checkForInf(neural.data, 'neural', options);
                if ~passed
                    message = 'neural.data contains Inf values';
                end
            else
                passed = true;
            end
            
        % Compatibility checks
        case 'fsMatch'
            if ~isempty(stim) && ~isempty(neural) && ...
               isfield(stim, 'fs') && isfield(neural, 'fs')
                passed = stim.fs == neural.fs;
                if ~passed
                    message = sprintf('Sampling frequency mismatch: stim=%g Hz, neural=%g Hz', ...
                                    stim.fs, neural.fs);
                end
            else
                passed = true;
            end
            
        case 'nTrialsMatch'
            if ~isempty(stim) && ~isempty(neural) && ...
               isfield(stim, 'data') && isfield(neural, 'data')
                stimNTrials = size(stim.data, 2);
                neuralNTrials = length(neural.data);
                passed = stimNTrials == neuralNTrials;
                if ~passed
                    message = sprintf('Trial count mismatch: stim=%d, neural=%d', ...
                                    stimNTrials, neuralNTrials);
                end
            else
                passed = true;
            end
            
        case 'trialLengthsCompat'
            if ~isempty(stim) && ~isempty(neural) && ...
               isfield(stim, 'data') && isfield(neural, 'data')
                passed = checkTrialLengths(stim, neural, options);
                if ~passed
                    message = 'Trial length mismatches found';
                end
            else
                passed = true;
            end
            
        otherwise
            % Unknown check - treat as passed with warning
            passed = true;
            message = sprintf('Unknown check: %s', checkName);
            severity = 'warning';
    end
end

function passed = checkStimDimensions(stim, options)
% Check stimulus feature dimensions consistency
    
    passed = true;
    [nFeatures, nTrials] = size(stim.data);
    
    if options.checkAllTrials
        trialsToCheck = 1:nTrials;
    else
        trialsToCheck = 1:min(3, nTrials);
    end
    
    try
        for f = 1:nFeatures
            if ~isempty(stim.data{f,1})
                featDim = size(stim.data{f,1}, 2);
                for t = trialsToCheck
                    if size(stim.data{f,t}, 2) ~= featDim
                        passed = false;
                        return;
                    end
                end
            end
        end
    catch
        passed = false;
    end
end

function passed = checkNeuralChannels(neural, options)
% Check neural channel consistency
    
    passed = true;
    
    if options.checkAllTrials
        trialsToCheck = 1:length(neural.data);
    else
        trialsToCheck = 1:min(5, length(neural.data));
    end
    
    try
        nChans = size(neural.data{1}, 2);
        for t = trialsToCheck(2:end) % Start from 2
            if size(neural.data{t}, 2) ~= nChans
                passed = false;
                return;
            end
        end
    catch
        passed = false;
    end
end

function passed = checkForNaN(dataCell, dataType, options)
% Check for NaN values in data
    
    passed = true;
    
    try
        if strcmp(dataType, 'stimulus')
            % For stimulus data (cell array with features x trials)
            nSamples = min(options.maxSamples, numel(dataCell));
            for idx = 1:nSamples
                if any(isnan(dataCell{idx}(:)))
                    passed = false;
                    return;
                end
            end
        else
            % For neural data (cell array with trials)
            nSamples = min(options.maxSamples, length(dataCell));
            for t = 1:nSamples
                if any(isnan(dataCell{t}(:)))
                    passed = false;
                    return;
                end
            end
        end
    catch
        passed = true; % Skip if can't check
    end
end

function passed = checkForInf(dataCell, dataType, options)
% Check for Inf values in data
    
    passed = true;
    
    try
        if strcmp(dataType, 'stimulus')
            % For stimulus data
            nSamples = min(options.maxSamples, numel(dataCell));
            for idx = 1:nSamples
                if any(isinf(dataCell{idx}(:)))
                    passed = false;
                    return;
                end
            end
        else
            % For neural data
            nSamples = min(options.maxSamples, length(dataCell));
            for t = 1:nSamples
                if any(isinf(dataCell{t}(:)))
                    passed = false;
                    return;
                end
            end
        end
    catch
        passed = true; % Skip if can't check
    end
end

function passed = checkTrialLengths(stim, neural, options)
% Check trial length compatibility
    
    passed = true;
    nTrials = min(size(stim.data, 2), length(neural.data));
    
    if options.checkAllTrials
        trialsToCheck = 1:nTrials;
    else
        trialsToCheck = 1:min(5, nTrials);
    end
    
    try
        mismatches = 0;
        for t = trialsToCheck
            if ~isempty(stim.data{1,t}) && ~isempty(neural.data{t})
                stimLen = size(stim.data{1,t}, 1);
                neuralLen = size(neural.data{t}, 1);
                % Check with tolerance
                if abs(stimLen - neuralLen) > max(stimLen, neuralLen) * options.tolerance
                    mismatches = mismatches + 1;
                    if mismatches >= 2 % Allow one mismatch as warning
                        passed = false;
                        return;
                    end
                end
            end
        end
    catch
        passed = false;
    end
end

function suggestion = getSuggestion(checkName, message)
% Get suggestion for fixing failed check
    
    switch checkName
        case {'stimIsStruct', 'neuralIsStruct'}
            suggestion = 'Ensure the file contains a valid MATLAB structure.';
            
        case {'stimHasData', 'neuralHasData'}
            suggestion = 'The data field is required. Add it as a cell array.';
            
        case {'stimHasFs', 'neuralHasFs'}
            suggestion = 'The fs (sampling frequency) field is required. Add it as a positive number.';
            
        case {'stimHasNames', 'neuralHasChanlocs'}
            suggestion = 'This field is marked as required in options. Add it or set requireStimNames/requireChanlocs to false.';
            
        case {'stimDataIsCell', 'neuralDataIsCell'}
            suggestion = 'Convert data to cell array format: use {} instead of [].';
            
        case {'stimDataNotEmpty', 'neuralDataNotEmpty'}
            suggestion = 'Ensure data contains at least one trial.';
            
        case {'stimDataNumeric', 'neuralDataNumeric'}
            suggestion = 'Ensure all data values are numeric (not strings or cells).';
            
        case {'stimFsPositive', 'neuralFsPositive'}
            suggestion = 'Sampling frequency must be > 0.';
            
        case 'fsMatch'
            suggestion = 'Ensure both files have the same sampling frequency or resample one dataset.';
            
        case 'nTrialsMatch'
            suggestion = 'Check trial counts and remove mismatched trials if necessary.';
            
        case 'trialLengthsCompat'
            suggestion = 'Consider trimming or padding trials to equal lengths.';
            
        otherwise
            suggestion = 'Review the CND format specification for requirements.';
    end
    
    % Add specific suggestions based on message
    if contains(message, 'NaN')
        suggestion = [suggestion ' Remove or replace NaN values with appropriate values.'];
    elseif contains(message, 'Inf')
        suggestion = [suggestion ' Remove or replace Inf values with appropriate values.'];
    end
end

function summaryTable = createSummaryTable(checkResults)
% Create summary table from results
    
    categoryNames = fieldnames(checkResults.categories);
    nCategories = length(categoryNames);
    
    % Create table data
    catNames = cell(nCategories, 1);
    passedCounts = zeros(nCategories, 1);
    failedCounts = zeros(nCategories, 1);
    warningCounts = zeros(nCategories, 1);
    totalCounts = zeros(nCategories, 1);
    
    for i = 1:nCategories
        catName = categoryNames{i};
        catStruct = checkResults.categories.(catName);
        
        % Convert from valid name back to display name
        displayName = strrep(catName, '_', ' ');
        displayName(1) = upper(displayName(1));
        
        catNames{i} = displayName;
        passedCounts(i) = catStruct.passed;
        failedCounts(i) = catStruct.failed;
        warningCounts(i) = length(catStruct.warnings);
        totalCounts(i) = catStruct.total;
    end
    
    % Add totals row
    catNames{end+1} = 'TOTAL';
    passedCounts(end+1) = checkResults.checksPassed;
    failedCounts(end+1) = checkResults.checksTotal - checkResults.checksPassed;
    warningCounts(end+1) = sum(cellfun(@(x) length(x.warnings),struct2cell(checkResults.categories)));
    totalCounts(end+1) = checkResults.checksTotal;
    
    % Create table
    summaryTable = table(catNames, passedCounts, warningCounts, failedCounts, totalCounts, ...
        'VariableNames', {'Category', 'Passed', 'Warnings', 'Failed', 'Total'});
end

function checkResults = calculateValidationFlags(checkResults)
% Calculate final validation flags
    
    % Initialize flags
    stimFailed = false;
    neuralFailed = false;
    compatibilityFailed = false;
    
    % Get category names
    categoryNames = fieldnames(checkResults.categories);
    
    % Check each category for failures
    for i = 1:length(categoryNames)
        catName = categoryNames{i};
        catStruct = checkResults.categories.(catName);
        
        % Check if category has any failures
        if catStruct.failed > 0
            % Determine what type of checks are in this category
            switch catName
                case {'StructureIntegrity', 'RequiredFields', 'DataFormat', 'DataConsistency'}
                    % These categories contain both stimulus and neural checks
                    % We'll check individual checks below
                case 'Compatibility'
                    compatibilityFailed = true;
            end
        end
    end
    
    % Check individual checks for stimulus/neural failures
    for i = 1:length(checkResults.allChecks)
        check = checkResults.allChecks{i};
        checkName = check{3};  % e.g., 'stimIsStruct', 'neuralIsStruct'
        result = check{4};
        
        if ~result.passed && strcmp(result.severity, 'error')
            if contains(checkName, 'stim')
                stimFailed = true;
            elseif contains(checkName, 'neural')
                neuralFailed = true;
            end
        end
    end
    
    % Add flags to results
    checkResults.stimValid = ~stimFailed;
    checkResults.neuralValid = ~neuralFailed;
    checkResults.compatibilityValid = ~compatibilityFailed;
    checkResults.passedAll = (checkResults.checksPassed == checkResults.checksTotal);
end
function displayResults(summaryTable, checkResults)
% Display results in verbose mode
    
    fprintf('\n=== VALIDATION SUMMARY ===\n\n');
    disp(summaryTable);
    
    % Show data statistics
    if ~isempty(fieldnames(checkResults.dataStats))
        fprintf('\n=== DATA STATISTICS ===\n\n');
        if isfield(checkResults.dataStats, 'stim')
            fprintf('Stimulus Data:\n');
            fprintf('  Sampling frequency: %.1f Hz\n', checkResults.dataStats.stim.fs);
            fprintf('  Features: %d\n', checkResults.dataStats.stim.nFeatures);
            fprintf('  Trials: %d\n', checkResults.dataStats.stim.nTrials);
            if isfield(checkResults.dataStats.stim, 'trialLength')
                fprintf('  Trial length: %d samples\n', checkResults.dataStats.stim.trialLength);
                fprintf('  Features per trial: %d\n', checkResults.dataStats.stim.nFeaturesPerTrial);
            end
            fprintf('\n');
        end
        
        if isfield(checkResults.dataStats, 'neural')
            fprintf('Neural Data:\n');
            fprintf('  Sampling frequency: %.1f Hz\n', checkResults.dataStats.neural.fs);
            fprintf('  Trials: %d\n', checkResults.dataStats.neural.nTrials);
            if isfield(checkResults.dataStats.neural, 'trialLength')
                fprintf('  Trial length: %d samples\n', checkResults.dataStats.neural.trialLength);
                fprintf('  Channels: %d\n', checkResults.dataStats.neural.nChannels);
            end
            fprintf('\n');
        end
    end
    
    % Show failed checks by category
    anyFailures = false;
    anyWarnings = false;
    categoryNames = fieldnames(checkResults.categories);
    
    for i = 1:length(categoryNames)
        catName = categoryNames{i};
        catStruct = checkResults.categories.(catName);
        if catStruct.failed > 0, anyFailures = true; end
        if ~isempty(catStruct.warnings), anyWarnings = true; end
    end
    
    if anyFailures
        fprintf('\n=== FAILED CHECKS (Errors) ===\n\n');
        for i = 1:length(categoryNames)
            catName = categoryNames{i};
            catStruct = checkResults.categories.(catName);
            
            if ~isempty(catStruct.failures)
                displayName = strrep(catName, '_', ' ');
                displayName(1) = upper(displayName(1));
                fprintf('%s:\n', displayName);
                for j = 1:length(catStruct.failures)
                    fprintf('  • %s\n', catStruct.failures{j});
                end
                fprintf('\n');
            end
        end
    end
    
    if anyWarnings
        fprintf('\n=== WARNINGS ===\n\n');
        for i = 1:length(categoryNames)
            catName = categoryNames{i};
            catStruct = checkResults.categories.(catName);
            
            if ~isempty(catStruct.warnings)
                displayName = strrep(catName, '_', ' ');
                displayName(1) = upper(displayName(1));
                fprintf('%s:\n', displayName);
                for j = 1:length(catStruct.warnings)
                    fprintf('  • %s\n', catStruct.warnings{j});
                end
                fprintf('\n');
            end
        end
    end
    
    if ~anyFailures && ~anyWarnings
        fprintf('\n✅ All checks passed!\n\n');
    elseif anyFailures
        fprintf('\n❌ Validation failed with errors\n\n');
    else
        fprintf('\n⚠️  Validation passed with warnings\n\n');
    end
    
    % Show final status
    fprintf('Overall Status:\n');
    fprintf('  Stimulus: %s\n', getStatusText(checkResults.stimValid));
    fprintf('  Neural: %s\n', getStatusText(checkResults.neuralValid));
    if ~isempty(checkResults.compatibilityValid)
        fprintf('  Compatibility: %s\n', getStatusText(checkResults.compatibilityValid));
    end
    fprintf('  All Checks: %s (%d/%d passed)\n\n', ...
            getStatusText(checkResults.passedAll), ...
            checkResults.checksPassed, checkResults.checksTotal);
end

function statusText = getStatusText(isValid)
% Get status text for display
    if isValid
        statusText = '✅ PASS';
    else
        statusText = '❌ FAIL';
    end
end