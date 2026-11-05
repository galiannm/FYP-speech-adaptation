%% Create a super subject 
% Output: datasets2/LalorNatSpeech/dataCND/dataSub100.mat

clear; 
clc;

sourceMode = 'preprocessed'; % 'raw' or 'preprocessed'
subjectNumber = 100;
rootDir = fullfile('datasets2', 'speechEmotionAdultEEG', 'dataCND');

assert(isfolder(rootDir), 'Data folder not found: %s', rootDir);
assert(isfile(fullfile(rootDir, 'dataStim.mat')), ...
    'Shared dataStim.mat is required to ensure same stimulus per trial.');
Sstim = load(fullfile(rootDir, 'dataStim.mat'), 'stim');
nStimTrials = size(Sstim.stim.data, 2);

switch sourceMode
    case 'raw'
        sourceDir = rootDir;
        inputStem = 'dataSub';
        outputStem = 'dataSub';
    case 'preprocessed'
        sourceDir = fullfile(rootDir, '1-8Hz');
        inputStem = 'pre_dataSub';
        outputStem = 'pre_dataSub';
    otherwise
        error('Unknown sourceMode: %s', sourceMode);
end

assert(isfolder(sourceDir), 'Source folder missing: %s', sourceDir);
outFile = fullfile(sourceDir, sprintf('%s%d.mat', outputStem, subjectNumber));
assert(~isfile(outFile), 'Output already exists; not overwriting: %s', outFile);

allFiles = dir(fullfile(sourceDir, [inputStem '*.mat']));
subjectFiles = {};
subjectIds = [];
for k = 1:numel(allFiles)
    % extract subject number
    tok = regexp(allFiles(k).name, ['^' regexptranslate('escape', inputStem) '(\d+)\.mat$'], 'tokens', 'once');
    
    if isempty(tok), continue; end

    sid = str2double(tok{1}); % convert to num
    if sid == subjectNumber, continue; end

    % saving filename and subject Id
    subjectFiles{end+1} = fullfile(sourceDir, allFiles(k).name);  
    subjectIds(end+1) = sid;  
end
assert(~isempty(subjectFiles), 'No input subjects found in %s', sourceDir);
[subjectIds, ord] = sort(subjectIds);
subjectFiles = subjectFiles(ord);


sums = cell(1, nStimTrials); % running sum of EEG recordings for each trial
sampleCounts = cell(1, nStimTrials); % how many subjects contrinuted to each time sample
contributorsPerTrial = zeros(1, nStimTrials); % how many subjects contrinuted to each trial
minSubjectsPerSample = 2; % Discard trailing samples supported by <2 subjects.
acceptedIds = []; % which subject was included
reference = struct(); % store EEG config

for k = 1:numel(subjectFiles)
    sid = subjectIds(k);
    S = load(subjectFiles{k}, 'eeg');
    if ~isfield(S, 'eeg') || ~isstruct(S.eeg)
        warning('Subject %d has no valid eeg struct: skipped', sid);
        continue
    end
    subj = S.eeg;

    if ~isfield(subj, 'data') || ~iscell(subj.data) || ...
       ~isfield(subj, 'fs') || ~isfield(subj, 'chanlocs') || ...
       isempty(subj.chanlocs)
        warning('Subject %d is missing data/fs/chanlocs: skipped', sid);
        continue
    end

    labels = strtrim(string({subj.chanlocs.labels}));
    if numel(unique(labels)) ~= numel(labels)
        warning('Subject %d contains duplicate channel labels: skipped', sid);
        continue
    end
    
    % use the first valid subject as reference
    % to keep subsequent subjects consistent
    if isempty(fieldnames(reference))
        reference.fs = double(subj.fs);
        reference.labels = labels;
        reference.chanlocs = subj.chanlocs;
        reference.nChannels = numel(labels);
        if isfield(subj, 'paddingStartSample')
            reference.paddingStartSample = subj.paddingStartSample;
        end
    end

    if double(subj.fs) ~= reference.fs
        warning('Subject %d has mismatched sampling frequency: skipped', sid);
        continue
    end

    [sameLabels, orderInSubject] = ismember(reference.labels, labels);
    if ~all(sameLabels) || numel(labels) ~= reference.nChannels
        warning('Subject %d has incompatible electrode labels: skipped', sid);
        continue
    end

    hasPadding = isfield(subj, 'paddingStartSample');
    referenceHasPadding = isfield(reference, 'paddingStartSample');
    if hasPadding ~= referenceHasPadding || ...
       (hasPadding && ~isequal(subj.paddingStartSample, reference.paddingStartSample))
        warning('Subject %d has incompatible stimulus-onset padding: skipped', sid);
        continue
    end


    nData = numel(subj.data); % num trials

    % check if all trials are present
    if nData == nStimTrials 
        trialMap = 1:nStimTrials;
    
    % not all trials are present
    elseif isfield(subj, 'nonemptyTrials')
        nt = subj.nonemptyTrials(:)';
        if islogical(nt) && numel(nt) == nStimTrials && nnz(nt) == nData
            trialMap = find(nt);

        elseif isnumeric(nt) && numel(nt) == nData && ...
               all(isfinite(nt)) && all(nt == round(nt)) && ...
               all(nt >= 1) && all(nt <= nStimTrials) && ...
               numel(unique(nt)) == nData
            trialMap = nt;

        elseif isnumeric(nt) && numel(nt) == nStimTrials && ...
               all(ismember(nt, [0, 1])) && nnz(nt) == nData
            trialMap = find(nt);

        else
            warning('Subject %d has ambiguous trial mapping: skipped', sid);
            continue
        end
    else
        warning('Subject %d has %d of %d trials and no trial mapping: skipped', ...
            sid, nData, nStimTrials);
        continue
    end

    usedAnyTrial = false;
    for t = 1:nData
        stimTrial = trialMap(t); 
        X = subj.data{t}; % extract EEG recording
        if isempty(X), continue; end
        if ~isnumeric(X) || ndims(X) ~= 2 || size(X, 2) ~= reference.nChannels
            warning('Subject %d, trial %d: wrong EEG dimensions; omitted', sid, stimTrial);
            continue
        end
        if ~all(isfinite(X(:)))
            warning('Subject %d, trial %d: nonfinite EEG; omitted', sid, stimTrial);
            continue
        end

        X = double(X(:, orderInSubject));
        newLen = size(X, 1);

        if contributorsPerTrial(stimTrial) == 0
            sums{stimTrial} = zeros(size(X));
            sampleCounts{stimTrial} = zeros(newLen, 1);
        else
            oldLen = size(sums{stimTrial}, 1);
            if newLen > oldLen
                sums{stimTrial}(oldLen+1:newLen, :) = 0;
                sampleCounts{stimTrial}(oldLen+1:newLen, 1) = 0;
            end
        end
        
        % Add EEG recordings
        sums{stimTrial}(1:newLen, :) = sums{stimTrial}(1:newLen, :) + X;
        sampleCounts{stimTrial}(1:newLen, 1) = sampleCounts{stimTrial}(1:newLen, 1) + 1;
        contributorsPerTrial(stimTrial) = contributorsPerTrial(stimTrial) + 1;
        usedAnyTrial = true;
    end
    if usedAnyTrial
        acceptedIds(end+1) = sid; 
        fprintf('Included subject %d\n', sid);
    end
    clear S subj
end

if numel(acceptedIds) < 2
    error('Only %d subjects contributed. Need >=2 to compute a group average.', ...
        numel(acceptedIds));
end
if any(contributorsPerTrial == 0)
    error(['No EEG is available for stimulus trials: %s. ' ...
        'Review EEG-to-stimulus trial mapping before saving.'], ...
        mat2str(find(contributorsPerTrial == 0)));
end

% Create a clean CND-like EEG struct, without copying subject-specific
% reference/extChan/presentation-order metadata from first participant.
eeg = struct();
eeg.dataType = 'EEG';
eeg.deviceName = 'GrandAverageVirtualSubject';
eeg.fs = reference.fs;
eeg.chanlocs = reference.chanlocs;
eeg.extChan = {};
eeg.data = cell(1, nStimTrials);

for t = 1:nStimTrials
    lastValidSample = find(sampleCounts{t} >= minSubjectsPerSample, 1, 'last');
    if isempty(lastValidSample)
        error('Stimulus trial %d has no samples with %d contributors', t, minSubjectsPerSample);
    end

    % compute the average
    eeg.data{t} = sums{t}(1:lastValidSample,:) ./ sampleCounts{t}(1:lastValidSample);
end
if isfield(reference, 'paddingStartSample')
    eeg.paddingStartSample = reference.paddingStartSample;
end

eeg.nonemptyTrials = 1:nStimTrials;
eeg.sourceSubjectIds = acceptedIds;
eeg.subjectCountsPerTrial = contributorsPerTrial;
eeg.minSubjectsPerSample = minSubjectsPerSample;
eeg.averagingMethod = 'Unweighted samplewise mean of same stimulus-aligned EEG trials';

fprintf('\nContributors per trial: %s\n', mat2str(contributorsPerTrial));
fprintf('Accepted subjects: %s\n', mat2str(acceptedIds));
fprintf('EEG fs: %g Hz | channels: %d | trials: %d\n', ...
    eeg.fs, numel(eeg.chanlocs), numel(eeg.data));

save(outFile, 'eeg', '-v7.3');
fprintf('\nSaved super-subject to:\n%s\n', outFile);

% Quick numerical sanity check: one super-subject trial shape
fprintf('First trial: %d samples x %d channels\n', size(eeg.data{1},1), size(eeg.data{1},2));
