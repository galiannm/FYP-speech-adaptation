
%% CNSP Dataset Feature Inventory
% Scans all dataset folders and extracts stimulus feature names.
%
% Expected structure:
% datasets2/
%   DiliBach/
%       dataCND/
%           dataStim1.mat
%           dataStim2.mat
%   AAD_KULeuven/
%       dataCND/
%           dataStim1.mat
%           dataStim2.mat
%           ...
%
% Output:
% 1. datasets_features.csv -> FolderName | Features
% 2. features_datasets.csv -> FeatureName | Folders
% 3. grouped_features_datasets.csv -> FeatureGroup | Folders

clear;
clc;

%% CONFIGs

rootDir = '/Users/manon/Documents/Uni/Trinity4/FYP/CNSP-resources/CNSP/datasets2';
outputDir = '/Users/manon/Documents/Uni/Trinity4/FYP/CNSP-resources/CNSP/myOutputs';

if ~isfolder(rootDir)
    error('Dataset folder does not exist: %s', rootDir);
end

% Discover dataset folders
folders = dir(rootDir);
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.', '..'}));

datasetNames = strings(0,1);
datasetFeatures = {};

fprintf('\n========== EXTRACT FEATURES FROM DATASETS ==========\n');

for d = 1:length(folders)

    datasetName = string(folders(d).name);

    % Locate the dataCND folder
    dataPath = fullfile(rootDir, datasetName, 'dataCND');

    if ~isfolder(dataPath)
        fprintf('[SKIP] %s: no dataCND folder\n', datasetName);
        continue;
    end

    % Find all stimulus files
    stimFiles = dir(fullfile(dataPath, 'dataStim*.mat'));

    if isempty(stimFiles)
        fprintf('[SKIP] %s: no stimulus files\n', datasetName);
        continue;
    end

    fprintf('\nDataset: %s\n', datasetName);

    features = strings(0,1);

    for f = 1:length(stimFiles)

        filePath = fullfile(dataPath, stimFiles(f).name);

        try

            % Load only the stim variable
            S = load(filePath, 'stim');

            if ~isfield(S, 'stim') || ...
                    ~isfield(S.stim, 'names')

                warning('No stim.names in %s', filePath);
                continue;

            end

            % Extract names
            rawNames = S.stim.names;

            if ischar(rawNames)
                names = string(cellstr(rawNames));
            else
                names = string(rawNames(:));
            end
            
            names = names(:);

            % Clean names
            names = strtrim(names);
            names = names(~ismissing(names) & names ~= "");

            % Merge names without duplicates
            features = union(features, names);

        catch ME

            warning('Error reading %s: %s', ...
                filePath, ME.message);

        end

    end

    % Sort alphabetically
    features = sort(features);

    % Store dataset information
    datasetNames(end+1,1) = datasetName;
    datasetFeatures{end+1,1} = features;

    fprintf('  Features found: %d\n', numel(features));

end

nDatasets = numel(datasetNames);


%% DATASET -> LIST OF FEATURES

% One row per dataset.
% All features are joined into a single string.

featureLists = strings(nDatasets,1);

for d = 1:nDatasets

    features = datasetFeatures{d};
    featureLists(d) = strjoin(features, '; ');

end

% Create output table
datasetTable = table( ...
    datasetNames, ...
    featureLists, ...
    'VariableNames', {'FolderName', 'Features'} ...
);

% Sort alphabetically by folder name
datasetTable = sortrows(datasetTable, 'FolderName');

writetable(datasetTable, ...
    fullfile(outputDir, 'datasets_features.csv'));


%% CREATE FEATURE LIST

% Collect every unique original feature name

allFeatures = strings(0,1);

for d = 1:nDatasets

    allFeatures = union( ...
        allFeatures, datasetFeatures{d});

end

% Alphabetical sorting
[~, idx] = sort(lower(allFeatures));
allFeatures = allFeatures(idx);

nFeatures = numel(allFeatures);


%% FEATURE -> LIST OF DATASETS

% For each feature, find all datasets containing it

featureFolders = strings(nFeatures,1);

for i = 1:nFeatures

    currentFeature = allFeatures(i);

    presentIn = strings(0,1);

    for d = 1:nDatasets

        % Does this dataset contain the current feature?
        if ismember(currentFeature, datasetFeatures{d})

            presentIn(end+1,1) = datasetNames(d);

        end

    end

    % Store list of folders
    featureFolders(i) = strjoin(sort(presentIn), '; ');

end

% Create table
featureTable = table( ...
    allFeatures, ...
    featureFolders, ...
    'VariableNames', {'FeatureName', 'Folders'} ...
);

% Alphabetical sorting
featureTable = sortrows(featureTable, 'FeatureName');

% Export
writetable(featureTable, ...
    fullfile(outputDir, 'features_datasets.csv'));


%% GROUP SIMILAR FEATURE NAMES

% Map each original feature to a standardised name.
%
% IMPORTANT:
% These are proposed naming groups, not necessarily
% scientifically identical features.
%
% Unrecognised names are preserved as they are.

featureGroups = strings(nFeatures,1);
reviewNeeded = false(nFeatures,1);

for i = 1:nFeatures

    original = allFeatures(i);

    % Normalise spelling for matching
    name = lower(strtrim(original));

    name = regexprep(name, '[_\-]', ' ');
    name = regexprep(name, '\s+', ' ');

    % ========================================
    % ACOUSTIC FEATURES
    % ========================================

    if ismember(name, [ ...
            "env", ...
            "envelope", ...
            "envelopes", ...
            "acoustic envelope", ...
            "speech envelope", ...
            "speech envelope vectors" ...
            ])

        featureGroups(i) = "Envelope";

    elseif ismember(name, [ ...
            "env'", ...
            "envder", ...
            "env del", ...
            "envelope derivative" ...
            ])

        featureGroups(i) = "Envelope derivative";

    elseif ismember(name, [ ...
            "spectrogram", ...
            "sgram" ...
            ])

        featureGroups(i) = "Spectrogram";

    % ========================================
    % SPEECH / LINGUISTIC FEATURES
    % ========================================

    elseif ismember(name, [ ...
            "word onset", ...
            "word onsets", ...
            "wordonsets", ...
            "word onset vectors", ...
            "wo" ...
            ])

        featureGroups(i) = "Word onset";

    elseif ismember(name, [ ...
            "phoneme onset", ...
            "phoneme onsets", ...
            "phoneonsets" ...
            ])

        featureGroups(i) = "Phoneme onset";

    elseif name == "pho"

        featureGroups(i) = "Phoneme onset";
        reviewNeeded(i) = true;

    elseif ismember(name, [ ...
            "phonetic features", ...
            "phonetic features (step)", ...
            "phonfeatures" ...
            ])

        featureGroups(i) = "Phonetic features";

    elseif name == "surprisal"

        featureGroups(i) = "Surprisal";

    elseif name == "entropy"

        featureGroups(i) = "Entropy";

    elseif ismember(name, [ ...
            "semantic dissimilarity", ...
            "diss" ...
            ])

        featureGroups(i) = "Semantic dissimilarity";

        if name == "diss"
            reviewNeeded(i) = true;
        end

    % ========================================
    % MUSIC FEATURES
    % ========================================

    elseif ismember(name, [ ...
            "note onset", ...
            "note onsets" ...
            ])

        featureGroups(i) = "Note onset";

    elseif name == "note onsets (metronome excluded)"

        featureGroups(i) = "Note onset (metronome excluded)";

    elseif name == "note pitch vectors"

        featureGroups(i) = "Note pitch";

    elseif ismember(name, ["s pitch", "h pitch"])

        featureGroups(i) = "Note pitch (" + original + ")";
        reviewNeeded(i) = true;

    elseif ismember(name, ["s dur", "h dur"])

        featureGroups(i) = "Note duration (" + original + ")";
        reviewNeeded(i) = true;

    elseif name == "expectation vector (metronome excluded)"

        featureGroups(i) = ...
            "Musical expectation (metronome excluded)";

    % ========================================
    % UNKNOWN / DATASET-SPECIFIC FEATURES
    % ========================================

    else

        % Keep original name if no matching rule exists
        featureGroups(i) = original;

        % Flag internally for future manual review
        reviewNeeded(i) = true;

    end

end


%% GROUPED FEATURE -> LIST OF DATASETS

% Get all unique standardised groups
uniqueGroups = unique(featureGroups);

nGroups = numel(uniqueGroups);

groupFolders = strings(nGroups,1);

for g = 1:nGroups

    currentGroup = uniqueGroups(g);

    % Find original feature names belonging to this group
    groupFeatures = allFeatures( ...
        featureGroups == currentGroup);

    presentIn = strings(0,1);

    for d = 1:nDatasets

        % Find whether ANY feature belonging to this group
        % is present in the current dataset.

        if any(ismember( ...
                groupFeatures, datasetFeatures{d}))

            presentIn(end+1,1) = datasetNames(d);

        end

    end

    % Store sorted list of unique folders
    groupFolders(g) = strjoin( ...
        sort(unique(presentIn)), '; ');

end

% Create grouped feature table
groupTable = table( ...
    uniqueGroups, ...
    groupFolders, ...
    'VariableNames', {'FeatureGroup', 'Folders'} ...
);

% Alphabetical sorting
groupTable = sortrows(groupTable, 'FeatureGroup');

% Export
writetable(groupTable, ...
    fullfile(outputDir, 'grouped_features_datasets.csv'));


%% FINAL SUMMARY

fprintf('\n========================================\n');
fprintf('          SCANNING COMPLETE!\n');
fprintf('========================================\n');

fprintf('\nDatasets scanned: %d', nDatasets);
fprintf('\nOriginal feature names: %d', nFeatures);
fprintf('\nStandardised feature groups: %d', nGroups);
fprintf('\nFeatures flagged for review: %d\n', ...
    sum(reviewNeeded));

fprintf('\nGenerated files:\n');

fprintf('1. datasets_features.csv\n');
fprintf('2. features_datasets.csv\n');
fprintf('3. grouped_features_datasets.csv\n');

fprintf('\nSaved in:\n%s\n', outputDir);


%% SANITY CHECK — INSPECT A SPECIFIC DATASET

fprintf('\n========================================\n');
fprintf('         DATASET SANITY CHECK\n');
fprintf('========================================\n');

dataset = input('Enter dataset folder name: ', 's');

dataPath = fullfile(rootDir, dataset, 'dataCND');
stimFiles = dir(fullfile(dataPath, 'dataStim*.mat'));

fprintf('\n=== Features in %s ===\n', dataset);

for f = 1:numel(stimFiles)

    S = load(fullfile(dataPath, stimFiles(f).name), 'stim');

    fprintf('\n%s:\n', stimFiles(f).name);
    disp(string(S.stim.names(:)));

end