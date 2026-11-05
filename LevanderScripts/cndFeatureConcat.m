function [stimFeature, feaGroupIdxs] = cndFeatureConcat(stim, stimIdxs)
% cndFeatureConcat Concatenate stimulus features (optionally grouped).
%
% Inputs
%   stim     : stimulus struct with stim.data{feature,trial} and stim.names
%   stimIdxs : cell array of feature-index vectors defining feature groups
%
% Outputs
%   stimFeature  : stim struct with stimFeature.data{trial} = concatenated features
%   feaGroupIdxs : 1 x nCols group index for each output feature column
%
% Example
%   [s1,g1] = cndFeatureConcat(stim,{[1],[2],[3]});
%   [s2,g2] = cndFeatureConcat(stim,{[1,2],[3]});
%   % s1.data and s2.data are identical; g1 vs g2 reflect grouping

    if ~isstruct(stim) || ~isfield(stim, 'data') || ~iscell(stim.data)
        error('stim must be a struct with stim.data as a cell array.');
    end
    if ~iscell(stimIdxs) || isempty(stimIdxs)
        error('stimIdxs must be a non-empty cell array of index vectors.');
    end

    nFeats = size(stim.data, 1);
    nTrials = size(stim.data, 2);

    stimFeature = stim;
    stimFeature.data = cell(1, nTrials);

    feaGroupIdxs = [];

    for iGroup = 1:numel(stimIdxs)
        idxs = stimIdxs{iGroup};
        if isempty(idxs) || ~isnumeric(idxs) || any(idxs < 1) || any(idxs > nFeats) || any(mod(idxs,1) ~= 0)
            error('Each stimIdxs{%d} must contain valid feature indices.', iGroup);
        end

        for iLocal = 1:numel(idxs)
            feaIdx = idxs(iLocal);

            nCols = size(stim.data{feaIdx, 1}, 2);
            for tr = 1:nTrials
                if size(stim.data{feaIdx, tr}, 2) ~= nCols
                    error('Feature %d has inconsistent dimensionality across trials.', feaIdx);
                end
            end
            feaGroupIdxs = [feaGroupIdxs, iGroup * ones(1, nCols)]; %#ok<AGROW>

            for tr = 1:nTrials
                X = stim.data{feaIdx, tr};
                if isempty(stimFeature.data{tr})
                    stimFeature.data{tr} = X;
                else
                    minLen = min(size(stimFeature.data{tr}, 1), size(X, 1));
                    stimFeature.data{tr} = [stimFeature.data{tr}(1:minLen, :), X(1:minLen, :)];
                end
            end
        end
    end

    if isfield(stimFeature, 'names') && ~isempty(stimFeature.names)
        stimFeature.names = {'Combined feature'};
    end














