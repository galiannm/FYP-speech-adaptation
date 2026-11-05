function avgModel = mTRFmodelAvg(models, param)
%MTRFMODELAVG Calculate average TRF model across subjects (BACKWARDS COMPATIBLE)
%
% Purpose:
%   Average TRF weights across multiple subjects. Supports both old and new syntax.
%
% Inputs (BACKWARDS COMPATIBLE):
%   models - struct array of TRF models (one per subject)
%   param  - If scalar 0/1: normalization flag (old syntax - what plotTRFweightsSEM uses)
%            If vector: feature indices to average (new syntax)
%            If empty: default to no normalization and all features
%
% Outputs:
%   avgModel - averaged model structure
%
% Examples (ALL WORK):
%   % OLD SYNTAX (what plotTRFweightsSEM uses):
%   avgModel = mTRFmodelAvg(modelAll, 1);  % Normalize weights
%   avgModel = mTRFmodelAvg(modelAll, 0);  % Don't normalize
%
%   % NEW SYNTAX (feature selection):
%   avgModel = mTRFmodelAvg(modelAll, [1 3]);  % Average features 1 and 3
%   avgModel = mTRFmodelAvg(modelAll, 1:5);    % Average first 5 features
%
%   % DEFAULT: no normalization, all features
%   avgModel = mTRFmodelAvg(modelAll);

    % ============================================
    % BACKWARDS COMPATIBILITY LOGIC
    % ============================================
    
    % Handle missing parameter
    if nargin < 2
        param = [];
    end
    
    % DETERMINE WHICH SYNTAX IS BEING USED
    % Old syntax: param is 0 or 1 (normalization flag) - used by plotTRFweightsSEM
    % New syntax: param is vector of feature indices
    
    isOldSyntax = false;
    normFlag = 0;
    featureIdx = [];
    
    if isempty(param)
        % Default: use old syntax with no normalization
        isOldSyntax = true;
        normFlag = 0;
    elseif isscalar(param) && (param == 0 || param == 1)
        % OLD SYNTAX: normalization flag (0 or 1)
        % This is what plotTRFweightsSEM calls with!
        isOldSyntax = true;
        normFlag = param;
    else
        % NEW SYNTAX: feature indices
        isOldSyntax = false;
        featureIdx = param;
        
        % Validate feature indices
        if ~isnumeric(featureIdx) || ~isvector(featureIdx)
            error('Feature indices must be a numeric vector');
        end
    end
    
    nSubs = length(models);
    
    % ============================================
    % OLD SYNTAX PATH (what plotTRFweightsSEM expects)
    % ============================================
    if isOldSyntax
        % Get weights from first subject to determine dimensions
        w1 = models(1).w;
        
        % Handle 2D or 3D weight matrices
        if ndims(w1) == 3
            % 3D: [features × time × channels] - need to reshape for old syntax
            [nFeatures, nTimeLags, nChannels] = size(w1);
            
            % Reshape to 2D: [channels × time] for each subject
            wAvg = zeros(nChannels, nTimeLags);
            
            for sub = 1:nSubs
                % Extract weights and reshape: mean across features
                w = mean(models(sub).w, 1);  % Mean across features -> [1 × time × channels]
                w = squeeze(w);              % -> [time × channels]
                w = w';                      % -> [channels × time]
                
                if normFlag
                    w = w / std(w(:));
                end
                
                if sub == 1
                    wAvg = w;
                else
                    wAvg = wAvg + w;
                end
            end
            wAvg = wAvg / nSubs;
            
        else
            % 2D: [channels × time] - original old syntax
            wAvg = models(1).w;
            if normFlag
                wAvg = wAvg / std(wAvg(:));
            end
            
            for sub = 2:nSubs
                w = models(sub).w;
                if normFlag
                    w = w / std(w(:));
                end
                wAvg = wAvg + w;
            end
            wAvg = wAvg / nSubs;
        end
        
        % Create output structure (old format)
        avgModel = models(1);
        avgModel.w = wAvg;
        
    % ============================================
    % NEW SYNTAX PATH (feature selection)
    % ============================================
    else
        % Get dimensions from first subject
        w1 = models(1).w;
        maxFeatures = size(w1, 1);
        
        % Validate feature indices
        if any(featureIdx > maxFeatures) || any(featureIdx < 1)
            error('Invalid feature indices. Maximum feature index is %d', maxFeatures);
        end
        
        nFeatures = length(featureIdx);
        nTimeLags = size(w1, 2);
        
        % Handle 2D or 3D weight matrices
        if ndims(w1) == 3
            nChannels = size(w1, 3);
            w_acc = zeros(nFeatures, nTimeLags, nChannels);
            
            for sub = 1:nSubs
                w_acc = w_acc + models(sub).w(featureIdx, :, :);
            end
            w_avg = w_acc / nSubs;
        else
            % 2D weights (backwards compatibility)
            w_acc = zeros(nFeatures, nTimeLags);
            
            for sub = 1:nSubs
                w_acc = w_acc + models(sub).w(featureIdx, :);
            end
            w_avg = w_acc / nSubs;
        end
        
        % Create enhanced output structure (new format)
        avgModel = struct();
        avgModel.w = w_avg;
        avgModel.t = models(1).t;
        
        if isfield(models(1), 'fs')
            avgModel.fs = models(1).fs;
        end
        
        avgModel.nSubs = nSubs;
        avgModel.featureIdx = featureIdx;
        
        % Copy optional fields
        if isfield(models(1), 'dir')
            avgModel.dir = models(1).dir;
        end
        
        if isfield(models(1), 'lambda')
            lambda_vals = zeros(1, nSubs);
            for sub = 1:nSubs
                if isfield(models(sub), 'lambda')
                    lambda_vals(sub) = models(sub).lambda;
                end
            end
            avgModel.lambda = mean(lambda_vals);
        end
        
        % Copy feature names if present
        if isfield(models(1), 'featureNames')
            avgModel.featureNames = models(1).featureNames(featureIdx);
        end
    end
end