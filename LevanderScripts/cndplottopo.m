function cndplottopology(data, chanlocs, varargin)
    % CNDPLOTTOPO: Interactive EEG Topographic Map Browser
    %
    % DESCRIPTION:
    %   Interactive visualization of EEG topographic maps with hover tooltips,
    %   click identification, and multi-frame navigation. Built on EEGLAB's
    %   topoplot function with enhanced interactivity.
    %
    % FEATURES:
    %   - Hover tooltips showing channel names and indices
    %   - Click identification to pin channel info in title
    %   - Automatic slider for multi-frame data browsing
    %   - Mouse scroll navigation through frames
    %   - Flexible electrode display options (dots, labels, both, numbers)
    %   - Customizable title, legend, and frame naming
    %
    % INPUTS:
    %   data      - EEG data matrix:
    %               Single frame:  [1 × N_channels] vector
    %               Multiple frames: [N_frames × N_channels] matrix
    %               Empty []: Display only electrode layout
    %   chanlocs  - EEGLAB channel locations structure
    %               Must contain .labels field (cell array of channel names)
    %               Can include .X/.Y or .theta/.radius coordinates
    %
    % OPTIONAL PARAMETER-VALUE PAIRS:
    %   'title'     - Plot title string (default: 'Topoplot')
    %   'legend'    - Colorbar label string (default: 'Value')
    %   'labels'    - Electrode display style:
    %                 'off'    : Show dots only (default)
    %                 'on'     : Show channel labels only
    %                 'both'   : Show both dots and labels
    %                 'numbers': Show channel indices instead of names
    %   'hover'     - Enable hover tooltips:
    %                 'on'  : Enable hover (default)
    %                 'off' : Disable hover
    %                 Note: Hover is automatically disabled for multi-frame data
    %   'frameName' - Label for frame counter (default: 'Frame')
    %                 e.g., 'Subject', 'Time', 'Condition'
    %
    % OUTPUTS:
    %   None directly returned. Creates an interactive figure with the
    %   specified EEG topography and interactive controls.
    %
    % USAGE EXAMPLES:
    %   % 1. Basic static plot
    %   cndplottopo(meanData, chanlocs, 'title', 'Mean Activity');
    %
    %   % 2. Interactive plot with hover
    %   cndplottopo(data, chanlocs, 'hover', 'on', 'labels', 'both');
    %
    %   % 3. Multi-subject browser
    %   cndplottopo(allSubjectsData, chanlocs, 'frameName', 'Subject');
    %
    %   % 4. Channel layout only
    %   cndplottopo([], chanlocs, 'labels', 'on', 'title', 'EEG Layout');
    %
    % INTERACTIVE FEATURES:
    %   - Single-frame data: Hover over electrodes to see channel info
    %   - Multi-frame data: Use slider or mouse scroll to browse frames
    %   - Click any electrode to display its name/index in the title
    %   - Color scale can be adjusted with caxis() after function call
    %
    % DATA ORIENTATION:
    %   Function automatically transposes data if needed to ensure format
    %   [N_frames × N_channels]. If data is [N_channels × N_frames], it will
    %   be transposed automatically.
    %
    % DEPENDENCIES:
    %   - EEGLAB's topoplot function
    %   - Standard MATLAB graphics
    %
    % LIMITATIONS:
    %   - Hover functionality disabled for multi-frame data to prevent
    %     conflicts with slider interaction
    %   - Maximum recommended channels: 128 for optimal display
    %
    % VERSION: 1.0
    % DATE: 03/02/2026
    % AUTHOR: Evander
    
    p = inputParser;
    addParameter(p, 'title', 'Topoplot', @ischar);
    addParameter(p, 'legend', 'Value', @ischar);
    addParameter(p, 'labels', 'off', @ischar); 
    addParameter(p, 'hover', 'on', @ischar);
    addParameter(p, 'frameName', 'Frame', @ischar); 
    parse(p, varargin{:});
    
    % Data orientation and state initialization
    if isempty(data)
        state.data = zeros(1, length(chanlocs));
        state.totalFrames = 1;
        state.isBlank = true;
        state.isMovie = false;
    else
        if size(data, 1) == length(chanlocs) && size(data, 2) ~= length(chanlocs)
            data = data'; 
        end
        state.data = data;
        state.totalFrames = size(data, 1);
        state.isBlank = false;
        state.isMovie = (size(data, 1) > 1);
    end
    
    state.chanlocs = chanlocs;
    state.params = p.Results;
    state.idx = 1;
    
    % Create figure with interactive controls
    fig = figure('Color', 'w', 'Name', 'CND Topo Browser', ...
                 'WindowScrollWheelFcn', @scrollCallback);
    
    % Hover only works for single-frame plots
    if strcmpi(state.params.hover, 'on') && ~state.isMovie
        set(fig, 'WindowButtonMotionFcn', @hoverCallback);
        state.hoverEnabled = true;
    else
        state.hoverEnabled = false;
    end
    
    % Add slider for multi-frame data
    if state.totalFrames > 1
        state.slider = uicontrol('Style', 'slider', 'Units', 'normalized', ...
            'Position', [0.2 0.03 0.6 0.04], 'Min', 1, 'Max', state.totalFrames, ...
            'Value', 1, 'SliderStep', [1/(state.totalFrames-1), 1/(state.totalFrames-1)], ...
            'Callback', @sliderCallback);
    end
    
    % Create hover tooltip if enabled
    if state.hoverEnabled
        state.hoverLabel = annotation('textbox', [0, 0, 0.3, 0.04], ...
            'String', '', 'Color', 'r', 'FontWeight', 'bold', 'BackgroundColor', 'w', ...
            'EdgeColor', 'k', 'Visible', 'off', 'Interpreter', 'none', ...
            'FitBoxToText', 'on', 'Margin', 5);
    else
        state.hoverLabel = [];
    end
        
    guidata(fig, state);
    renderPlot(fig);
end

function renderPlot(fig)
    state = guidata(fig);
    ax = findobj(fig, 'Type', 'axes');
    if ~isempty(ax), delete(ax); end
    
    % Map label options to topoplot electrode style
    requestedElec = lower(state.params.labels);
    switch requestedElec
        case 'off',     elecArg = 'on'; 
        case 'on',      elecArg = 'labels';
        case 'both',    elecArg = 'labelpoint';
        case 'numbers', elecArg = 'numpoint';
        otherwise,      elecArg = 'on';
    end
    
    if state.isBlank, plotStyle = 'blank'; else, plotStyle = 'both'; end
    
    currentData = double(state.data(state.idx, :));
    
    % Clear previous electrode tracking data
    state.plotX = [];
    state.plotY = [];
    state.plotLabels = {};
    state.plotChanIndices = [];
    
    % Create topoplot
    topoplot(currentData, state.chanlocs, ...
             'style', plotStyle, ...
             'electrodes', elecArg, ...
             'chaninfo', []);
    
    % Extract electrode positions for hover/click functionality
    ax = gca;
    children = get(ax, 'Children');
    
    state.electrodeHandles = [];
    
    for i = 1:length(children)
        obj = children(i);
        objType = get(obj, 'Type');
        
        if strcmp(objType, 'line')
            % Check for electrode dots (lines with markers)
            marker = get(obj, 'Marker');
            if ~isempty(marker) && ~strcmp(marker, 'none')
                xpos = get(obj, 'XData');
                ypos = get(obj, 'YData');
                
                if ~isempty(xpos)
                    numDots = length(xpos);
                    
                    for j = 1:numDots
                        state.plotX = [state.plotX, xpos(j)];
                        state.plotY = [state.plotY, ypos(j)];
                        
                        % Electrode dots are in same order as chanlocs
                        if j <= length(state.chanlocs)
                            state.plotLabels{end+1} = state.chanlocs(j).labels;
                            state.plotChanIndices(end+1) = j;
                        else
                            state.plotLabels{end+1} = sprintf('Ch%d', j);
                            state.plotChanIndices(end+1) = j;
                        end
                    end
                    
                    state.electrodeHandles = [state.electrodeHandles, obj];
                    set(obj, 'ButtonDownFcn', {@clickCallback, fig});
                end
            end
            
        elseif strcmp(objType, 'text')
            % Electrode labels
            pos = get(obj, 'Position');
            labelText = get(obj, 'String');
            
            state.plotX = [state.plotX, pos(1)];
            state.plotY = [state.plotY, pos(2)];
            state.plotLabels{end+1} = labelText;
            state.electrodeHandles = [state.electrodeHandles, obj];
            
            % Find matching channel
            allLabels = {state.chanlocs.labels};
            chanIdx = find(strcmp(allLabels, labelText));
            if ~isempty(chanIdx)
                state.plotChanIndices(end+1) = chanIdx;
            else
                numIdx = str2double(labelText);
                if ~isnan(numIdx) && numIdx <= length(state.chanlocs)
                    state.plotChanIndices(end+1) = numIdx;
                else
                    state.plotChanIndices(end+1) = length(state.plotChanIndices) + 1;
                end
            end
            
            set(obj, 'ButtonDownFcn', {@clickCallback, fig});
            
        elseif strcmp(objType, 'surface') || strcmp(objType, 'patch')
            % Remove click callback from heatmap surface
            try
                set(obj, 'ButtonDownFcn', '');
            catch
            end
        end
    end
    
    title(sprintf('%s (%s %d/%d)', state.params.title, state.params.frameName, ...
          state.idx, state.totalFrames));
      
    if ~state.isBlank
        cb = colorbar;
        ylabel(cb, state.params.legend);
    end
    
    guidata(fig, state);
end

function clickCallback(src, ~, fig)
    state = guidata(fig);
    
    objType = get(src, 'Type');
    
    if strcmp(objType, 'line')
        x = get(src, 'XData');
        y = get(src, 'YData');
        
        if length(x) > 1
            ax = gca;
            cp = get(ax, 'CurrentPoint');
            clickX = cp(1,1);
            clickY = cp(1,2);
            
            dist = sqrt((x - clickX).^2 + (y - clickY).^2);
            [~, minIdx] = min(dist);
            x = x(minIdx);
            y = y(minIdx);
        end
        
    elseif strcmp(objType, 'text')
        pos = get(src, 'Position');
        x = pos(1);
        y = pos(2);
    else
        return;
    end
    
    if ~isempty(state.plotX)
        x = x(1);
        y = y(1);
        
        plotX = state.plotX(:)';
        plotY = state.plotY(:)';
        
        dist = sqrt((plotX - x).^2 + (plotY - y).^2);
        [~, nearestIdx] = min(dist);
        
        if nearestIdx <= length(state.plotLabels)
            chanLabel = state.plotLabels{nearestIdx};
            if nearestIdx <= length(state.plotChanIndices)
                chanIdx = state.plotChanIndices(nearestIdx);
                title(['Clicked: ', chanLabel, ' (Index: ', num2str(chanIdx), ')'], 'Color', 'b');
            else
                title(['Clicked: ', chanLabel], 'Color', 'b');
            end
        end
    end
end

function hoverCallback(src, ~)
    state = guidata(src);
    
    if state.isMovie || ~state.hoverEnabled
        return;
    end
    
    if isempty(state.plotX) || isempty(state.plotLabels)
        return;
    end
    
    ax = gca; 
    cp = get(ax, 'CurrentPoint');
    mouseX = cp(1,1);
    mouseY = cp(1,2);
    
    plotX = state.plotX(:)';
    plotY = state.plotY(:)';
    
    dist = sqrt((plotX - mouseX).^2 + (plotY - mouseY).^2);
    [minDist, nearestIdx] = min(dist);
    
    if minDist < 0.06 && nearestIdx <= length(state.plotLabels)
        figPos = get(src, 'CurrentPoint'); 
        figSize = get(src, 'Position');
        
        chanLabel = state.plotLabels{nearestIdx};
        
        if nearestIdx <= length(state.plotChanIndices)
            chanIdx = state.plotChanIndices(nearestIdx);
            hoverText = sprintf(' Channel: %s | Index: %d ', chanLabel, chanIdx);
        else
            hoverText = sprintf(' Channel: %s ', chanLabel);
        end
        
        set(state.hoverLabel, ...
            'String', hoverText, ...
            'Position', [figPos(1)/figSize(3) + 0.02, figPos(2)/figSize(4) + 0.02, 0.3, 0.04], ...
            'Visible', 'on');
    else
        set(state.hoverLabel, 'Visible', 'off');
    end
end

function sliderCallback(src, ~)
    fig = ancestor(src, 'figure');
    state = guidata(fig);
    state.idx = round(get(src, 'Value'));
    guidata(fig, state);
    renderPlot(fig);
end

function scrollCallback(src, event)
    state = guidata(src);
    if state.totalFrames <= 1, return; end
    state.idx = max(1, min(state.totalFrames, state.idx + event.VerticalScrollCount));
    if isfield(state, 'slider'), set(state.slider, 'Value', state.idx); end
    guidata(src, state);
    renderPlot(src);
end






%% OPTION E: Visualise electrode layouts across datasets

dataRoot = 'datasets2';

dataFolders = {'AliceSpeech', 'FDSpeech_L1', 'LalorNatSpeech', 'childStoriesEEG_Sysoeva', ...
    'cocktailAttSwitch', 'cocktailAttSwitch', 'fNIRS_podcast', 'hyperscanningAmir_Listening', ...
    'podcastListeningEEG', 'ComedyShowEEG', 'speechEmotionAdultEEG', 'speechTrustEEG', };

for d = 1:length(dataFolders)

    datasetName = dataFolders{d};

    % Find the first subject file
    neuralFolder = fullfile(dataRoot, datasetName, 'dataCND');
    files = dir(fullfile(neuralFolder, 'dataSub*.mat'));

    searchPaths = {
        fullfile(neuralFolder, '1-8Hz', 'pre_dataSub*.mat')
        fullfile(neuralFolder, 'dataSub*.mat')
    };

    files = [];

    for p = 1:length(searchPaths)
        files = dir(searchPaths{p});

        if ~isempty(files)
            break;
        end
    end

    if isempty(files)
        fprintf('No files found for %s\n', datasetName);
        continue;
    end

    % Load electrode information
    filePath = fullfile(files(1).folder, files(1).name);
    S = load(filePath);

    if isfield(S, 'eeg')
        neural = S.eeg;
    elseif isfield(S, 'neural')
        neural = S.neural;
    elseif isfield(S, 'EEG')
        neural = S.EEG;
    else
        fprintf('No EEG structure found for %s\n', datasetName);
        continue;
    end

    if ~isfield(neural, 'chanlocs')
        fprintf('No chanlocs for %s\n', datasetName);
        continue;
    end

    chanlocs = neural.chanlocs;
    labels = {chanlocs.labels};

    % Convert electrode coordinates to double
    coordinateFields = {'theta', 'radius', ...
        'X', 'Y', 'Z', ...
        'sph_theta', 'sph_phi', 'sph_radius'};
    
    for c = 1:numel(chanlocs)
        for f = 1:numel(coordinateFields)
    
            field = coordinateFields{f};
    
            if isfield(chanlocs, field) && ...
                    ~isempty(chanlocs(c).(field)) && ...
                    isnumeric(chanlocs(c).(field))
    
                chanlocs(c).(field) = double(chanlocs(c).(field));
    
            end
        end
    end

    fprintf('\n%s: %d electrodes\n', datasetName, length(labels));
    fprintf('Loaded from: %s\n', filePath);
    disp(labels');

    % Display electrode layout
    cndplottopology([], chanlocs, ...
        'title', datasetName, ...
        'labels', 'both', ...
        'hover', 'on');

    fprintf('\n%s: %d electrodes\n', datasetName, length(labels));
    fprintf('Pz available: %d\n', any(strcmpi(labels, 'Pz')));
    fprintf('CPz available: %d\n', any(strcmpi(labels, 'CPz')));

end