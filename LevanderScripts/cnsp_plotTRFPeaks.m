function [peaks, fig] = cnsp_plotTRFPeaks(modelAll, chanlocs, varargin)
% CNSP_PLOTTRFPEAKS - Detect TRF peaks and create interactive topography viewer
%
% This function detects peaks and troughs in Temporal Response Function (TRF)
% data and creates an interactive visualization where clicking anywhere on the
% TRF line shows the corresponding scalp topography.
%
% INPUTS:
%   modelAll   : TRF models struct array from multiple subjects
%   chanlocs   : EEG channel locations structure (from EEGLAB)
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'dim2plot'      : Feature dimension to analyze (default: 1)
%   'timeRange'     : Time range in ms [min max] (default: [-100 600])
%   'nPeaks'        : Number of positive peaks to detect (default: 3)
%   'nTroughs'      : Number of negative peaks to detect (default: 3)
%   'displayTopoplot' : Show interactive figure (default: true)
%   'figTitle'      : Custom figure title (default: '')
%
% OUTPUTS:
%   peaks : Structure with peak/trough information
%   fig   : Figure handle (empty if displayTopoplot=false)

%% Parse input arguments
p = inputParser;
addRequired(p, 'modelAll');
addRequired(p, 'chanlocs');
addParameter(p, 'dim2plot', 1, @isnumeric);
addParameter(p, 'timeRange', [-100 600], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'nPeaks', 3, @(x) isempty(x) || (isscalar(x) && x>=0));
addParameter(p, 'nTroughs', 3, @(x) isempty(x) || (isscalar(x) && x>=0));
addParameter(p, 'displayTopoplot', true, @islogical);
addParameter(p, 'figTitle', '', @ischar);
parse(p, modelAll, chanlocs, varargin{:});

%% Prepare TRF data structure
TRFtopo = prepareTRFData(modelAll, chanlocs, p.Results.dim2plot, p.Results.timeRange);

%% Detect peaks and troughs
peaks = detectPeaksAndTroughs(TRFtopo, p.Results.nPeaks, p.Results.nTroughs);

%% Create interactive plot if requested
fig = [];
if p.Results.displayTopoplot
    fig = createInteractivePlot(TRFtopo, peaks, p.Results.figTitle);
end

end

function TRFtopo = prepareTRFData(modelAll, chanlocs, dim2plot, timeRange)
% PREPARETRFDATA - Convert model data to TRFtopo structure
%
% This function handles different data formats and creates a standardized
% TRFtopo structure with topographies for each time point.

fprintf('\n=== Preparing TRF Data ===\n');

% Convert cell array to struct array if needed
if iscell(modelAll)
    modelAll = [modelAll{:}];
end

% Determine data orientation and create averaged model
if ndims(modelAll(1).w) == 3
    % Format: [features, timepoints, channels]
    nFeatures = size(modelAll(1).w, 1);
    nTimePoints = size(modelAll(1).w, 2);
    nChannels = size(modelAll(1).w, 3);
    
    fprintf('Data format: [features=%d, timepoints=%d, channels=%d]\n', ...
        nFeatures, nTimePoints, nChannels);
    
    % Average across subjects
    avgW = zeros(nFeatures, nTimePoints, nChannels);
    for i = 1:length(modelAll)
        avgW = avgW + modelAll(i).w;
    end
    avgW = avgW / length(modelAll);
    
    % Use specified feature
    dim2plot = min(dim2plot, nFeatures);
    topoData = squeeze(avgW(dim2plot, :, :)); % [timepoints × channels]
    
elseif ndims(modelAll(1).w) == 2
    % Check orientation - assume [channels × timepoints] is common
    if size(modelAll(1).w, 1) == length(chanlocs)
        nChannels = size(modelAll(1).w, 1);
        nTimePoints = size(modelAll(1).w, 2);
        fprintf('Data format: [channels=%d × timepoints=%d]\n', nChannels, nTimePoints);
        
        % Average across subjects
        avgW = zeros(nChannels, nTimePoints);
        for i = 1:length(modelAll)
            avgW = avgW + modelAll(i).w;
        end
        avgW = avgW / length(modelAll);
        
        % Transpose to [timepoints × channels]
        topoData = avgW';
    else
        error('Cannot determine data orientation: number of channels mismatch');
    end
else
    error('Unexpected data dimensions: model.w must be 2D or 3D');
end

% Create TRFtopo structure
TRFtopo = struct();
TRFtopo.times = modelAll(1).t(:)';  % Time points as row vector
TRFtopo.nTimePoints = length(TRFtopo.times);
TRFtopo.nChannels = length(chanlocs);
TRFtopo.chanlocs = chanlocs;
TRFtopo.topographies = topoData; % [timepoints × channels]
TRFtopo.meanTRF = mean(TRFtopo.topographies, 2);
TRFtopo.dim2plot = dim2plot;

% Restrict to specified time range
timeIdx = TRFtopo.times >= timeRange(1) & TRFtopo.times <= timeRange(2);
if ~any(timeIdx)
    error('No time points found in specified range [%d, %d] ms', timeRange(1), timeRange(2));
end

TRFtopo.times = TRFtopo.times(timeIdx);
TRFtopo.topographies = TRFtopo.topographies(timeIdx, :);
TRFtopo.meanTRF = TRFtopo.meanTRF(timeIdx);
TRFtopo.nTimePoints = sum(timeIdx);

fprintf('Time range: %.1f to %.1f ms (%d points)\n', ...
    min(TRFtopo.times), max(TRFtopo.times), TRFtopo.nTimePoints);

% Verify data has variation
testIdx = round(TRFtopo.nTimePoints/2);
testTopo = TRFtopo.topographies(testIdx, :);
fprintf('Data verification at %.1f ms: range=[%.4f, %.4f], var=%.6f\n', ...
    TRFtopo.times(testIdx), min(testTopo), max(testTopo), var(testTopo));

end

function peaks = detectPeaksAndTroughs(TRFtopo, nPeaksRequested, nTroughsRequested)
% DETECTPEAKSANDTROUGHS - Find peaks and troughs in TRF data
%
% Uses first derivative method to detect local maxima and minima,
% then sorts by magnitude and selects top N.

fprintf('\n=== Detecting Peaks and Troughs ===\n');

w = TRFtopo.meanTRF(:);  % Ensure column vector
t = TRFtopo.times(:);

% First derivative method for peak detection
dw = diff(w);
signChange = diff(sign(dw));

% Find maxima (+ to -) and minima (- to +)
peakIdx = find(signChange < 0) + 1;
troughIdx = find(signChange > 0) + 1;

% Get values and times
peakTimesAll = t(peakIdx);
peakValsAll = w(peakIdx);
troughTimesAll = t(troughIdx);
troughValsAll = w(troughIdx);

% Sort by absolute magnitude
[~, peakOrder] = sort(abs(peakValsAll), 'descend');
[~, troughOrder] = sort(abs(troughValsAll), 'descend');

% Take top N as specified
nPeaksFound = min(nPeaksRequested, length(peakOrder));
nTroughsFound = min(nTroughsRequested, length(troughOrder));

if nPeaksFound > 0
    peakIdx_top = peakOrder(1:nPeaksFound);
else
    peakIdx_top = [];
end

if nTroughsFound > 0
    troughIdx_top = troughOrder(1:nTroughsFound);
else
    troughIdx_top = [];
end

% Sort by time for display (chronological order)
if ~isempty(peakIdx_top)
    [~, timeOrderPeaks] = sort(peakTimesAll(peakIdx_top));
    peakIdx_top = peakIdx_top(timeOrderPeaks);
    peakTimes = peakTimesAll(peakIdx_top);
    peakVals = peakValsAll(peakIdx_top);
else
    peakTimes = [];
    peakVals = [];
end

if ~isempty(troughIdx_top)
    [~, timeOrderTroughs] = sort(troughTimesAll(troughIdx_top));
    troughIdx_top = troughIdx_top(timeOrderTroughs);
    troughTimes = troughTimesAll(troughIdx_top);
    troughVals = troughValsAll(troughIdx_top);
else
    troughTimes = [];
    troughVals = [];
end

fprintf('Detected %d peaks and %d troughs\n', length(peakTimes), length(troughTimes));

% Extract topography for each detected event
allTimes = [peakTimes; troughTimes];
nTotal = length(allTimes);

if nTotal > 0
    topographyMatrix = zeros(nTotal, TRFtopo.nChannels);
    timepoints = zeros(1, nTotal);
    labels = cell(1, nTotal);
    
    for i = 1:nTotal
        [~, timeIdx] = min(abs(TRFtopo.times - allTimes(i)));
        topographyMatrix(i, :) = TRFtopo.topographies(timeIdx, :);
        timepoints(i) = TRFtopo.times(timeIdx);
        
        if i <= length(peakTimes)
            labels{i} = sprintf('P%d', i);
        else
            labels{i} = sprintf('N%d', i - length(peakTimes));
        end
    end
else
    topographyMatrix = [];
    timepoints = [];
    labels = {};
end

% Create output structure
peaks = struct();
peaks.positive.times = peakTimes;
peaks.positive.values = peakVals;
peaks.negative.times = troughTimes;
peaks.negative.values = troughVals;
peaks.topographies = topographyMatrix;
peaks.timepoints = timepoints;
peaks.labels = labels;
peaks.nPeaks = length(peakTimes);
peaks.nTroughs = length(troughTimes);
peaks.nTotal = nTotal;

% Display detected peaks/troughs
for i = 1:length(peakTimes)
    fprintf('  Peak %d: %.1f ms, value=%.3f\n', i, peakTimes(i), peakVals(i));
end
for i = 1:length(troughTimes)
    fprintf('  Trough %d: %.1f ms, value=%.3f\n', i, troughTimes(i), troughVals(i));
end

end

function fig = createInteractivePlot(TRFtopo, peaks, figTitle)
% CREATEINTERACTIVEPLOT - Create figure with interactive TRF topography viewer
%
% Creates a two-panel figure:
%   Left: TRF time series with detected peaks/troughs
%   Right: Topography map that updates when clicking on left plot
%
% Interactive features:
%   - Click anywhere on TRF line → shows topography at that time
%   - Mouse hover → shows vertical line with time indicator
%   - Click on detected peaks/troughs → shows their topographies

fprintf('\n=== Creating Interactive Figure ===\n');

% Create figure
fig = figure('Color', 'w', 'Position', [100, 100, 1400, 600], ...
    'Name', 'TRF Peak Analysis - Click anywhere on TRF line', ...
    'NumberTitle', 'off');

% Left subplot: TRF with peaks and troughs
ax1 = subplot(1, 2, 1);
hold on;

% Plot TRF line
hLine = plot(TRFtopo.times, TRFtopo.meanTRF, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'TRF');

% Plot peaks if they exist
if ~isempty(peaks.positive.times)
    hPeaks = plot(peaks.positive.times, peaks.positive.values, 'rv', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', 'Peaks');
end

% Plot troughs if they exist
if ~isempty(peaks.negative.times)
    hTroughs = plot(peaks.negative.times, peaks.negative.values, 'g^', ...
        'MarkerFaceColor', 'g', 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', 'Troughs');
end

% Add peak labels
for i = 1:length(peaks.positive.times)
    text(peaks.positive.times(i), peaks.positive.values(i), sprintf('P%d', i), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'r', ...
        'HandleVisibility', 'off');
end

% Add trough labels
for i = 1:length(peaks.negative.times)
    text(peaks.negative.times(i), peaks.negative.values(i), sprintf('N%d', i), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'g', ...
        'HandleVisibility', 'off');
end

% Formatting
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Mean TRF Weight', 'FontSize', 12);

if isempty(figTitle)
    title('Click anywhere on TRF line to see topography', ...
        'FontSize', 14, 'FontWeight', 'bold');
else
    title(figTitle, 'FontSize', 14, 'FontWeight', 'bold');
end

grid on;
xlim([TRFtopo.times(1), TRFtopo.times(end)]);

% Add zero lines
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, ...
    'HandleVisibility', 'off');
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, ...
    'HandleVisibility', 'off');

% Create legend - only show main elements
if ~isempty(peaks.positive.times) && ~isempty(peaks.negative.times)
    legend([hLine, hPeaks(1), hTroughs(1)], {'TRF', 'Peaks', 'Troughs'}, ...
        'Location', 'best', 'FontSize', 10);
elseif ~isempty(peaks.positive.times)
    legend([hLine, hPeaks(1)], {'TRF', 'Peaks'}, 'Location', 'best', 'FontSize', 10);
elseif ~isempty(peaks.negative.times)
    legend([hLine, hTroughs(1)], {'TRF', 'Troughs'}, 'Location', 'best', 'FontSize', 10);
else
    legend(hLine, {'TRF'}, 'Location', 'best', 'FontSize', 10);
end

% Right subplot: Topography placeholder
ax2 = subplot(1, 2, 2);
axis off;
text(0.5, 0.5, 'Click on left plot to see topography', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Color', [0.5 0.5 0.5]);

% Store data for callbacks
data.TRFtopo = TRFtopo;
data.peaks = peaks;
data.ax1 = ax1;
data.ax2 = ax2;
data.fig = fig;
data.currentHighlight = [];
data.vline = [];
data.vlineText = [];
guidata(fig, data);

% Set callbacks - make entire axes clickable
set(ax1, 'ButtonDownFcn', @trfClickCallback);
set(findobj(ax1, 'Type', 'line'), 'ButtonDownFcn', @trfClickCallback);
set(fig, 'WindowButtonMotionFcn', @trfMouseMoveCallback);
set(fig, 'WindowButtonUpFcn', @trfMouseUpCallback);

fprintf('Interactive figure ready:\n');
fprintf('  • Click ANYWHERE on TRF line → shows topography at that time\n');
fprintf('  • Hover mouse → shows time indicator line\n');
fprintf('  • Red triangles (P#) mark detected peaks\n');
fprintf('  • Green triangles (N#) mark detected troughs\n');

end

function trfClickCallback(src, ~)
% TRFCLICKCALLBACK - Handle clicks on TRF plot
%
% When user clicks on the TRF plot, this function:
%   1. Finds the closest time point to the click
%   2. Updates the right panel with topography at that time
%   3. Highlights the clicked point on the TRF plot

% Get figure data
data = guidata(gcf);
if isempty(data), return; end

% Get click position in axes coordinates
clickPoint = get(data.ax1, 'CurrentPoint');
clickTime = clickPoint(1,1);

% Find closest time point in TRF data
[~, timeIdx] = min(abs(data.TRFtopo.times - clickTime));
actualTime = data.TRFtopo.times(timeIdx);

% Get topography at this time
topoData = data.TRFtopo.topographies(timeIdx, :);

fprintf('Clicked at %.1f ms → showing topography at %.1f ms\n', ...
    clickTime, actualTime);
fprintf('  Topo stats: min=%.4f, max=%.4f, mean=%.4f\n', ...
    min(topoData), max(topoData), mean(topoData));

% Update topography display
updateTopographyPlot(data.ax2, topoData, data.TRFtopo.chanlocs, actualTime);

% Highlight the clicked point on TRF plot
axes(data.ax1);
hold on;

% Remove previous highlight if exists
if ~isempty(data.currentHighlight) && ishandle(data.currentHighlight)
    delete(data.currentHighlight);
end

% Add new highlight
data.currentHighlight = plot(actualTime, data.TRFtopo.meanTRF(timeIdx), 'ro', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r', ...
    'HandleVisibility', 'off');
hold off;

% Update stored data
guidata(gcf, data);
end

function trfMouseMoveCallback(~, ~)
% TRFMOUSEMOVECALLBACK - Show vertical line at mouse position
%
% When mouse moves over the TRF plot, this function:
%   1. Shows a vertical line at the current mouse x-position
%   2. Displays the time value at the top of the line

% Get figure data
data = guidata(gcf);
if isempty(data) || ~isfield(data, 'ax1'), return; end

% Get mouse position in axes coordinates
cp = get(data.ax1, 'CurrentPoint');
mouseX = cp(1,1);

% Check if mouse is within axes limits
xlims = get(data.ax1, 'XLim');
ylims = get(data.ax1, 'YLim');

if mouseX < xlims(1) || mouseX > xlims(2)
    % Mouse outside axes - hide the line
    if ~isempty(data.vline) && ishandle(data.vline)
        set(data.vline, 'Visible', 'off');
        set(data.vlineText, 'Visible', 'off');
    end
    return;
end

% Find closest time point
[~, timeIdx] = min(abs(data.TRFtopo.times - mouseX));
actualTime = data.TRFtopo.times(timeIdx);

% Update or create vertical line
axes(data.ax1);
hold on;

if isempty(data.vline) || ~ishandle(data.vline)
    % Create new vertical line
    data.vline = plot([actualTime actualTime], ylims, 'r--', ...
        'LineWidth', 1.5, 'Color', [0.8 0.2 0.2 0.7], ...
        'HandleVisibility', 'off');
    
    % Add time label
    data.vlineText = text(actualTime, ylims(2)*0.95, ...
        sprintf('%.1f ms', actualTime), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.8], ...
        'EdgeColor', [0.8 0.2 0.2], ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'HandleVisibility', 'off');
else
    % Update existing line
    set(data.vline, 'XData', [actualTime actualTime], 'Visible', 'on');
    set(data.vlineText, 'Position', [actualTime, ylims(2)*0.95], ...
        'String', sprintf('%.1f ms', actualTime), 'Visible', 'on');
end
hold off;

guidata(gcf, data);
end

function trfMouseUpCallback(~, ~)
% TRFMOUSEUPCALLBACK - Clean up after mouse release (placeholder)
%
% Currently empty but kept for future extensions (e.g., clicking and dragging)

% Optional: Add functionality here if needed
% For now, just ensure data is still valid
data = guidata(gcf);
if isempty(data), return; end

% Could add click-and-drag functionality here in the future
end

function updateTopographyPlot(ax, topoData, chanlocs, timepoint)
% UPDATETOPOGRAPHYPLOT - Update topography plot with fallback options
%
% This function attempts to display topography using multiple methods:
%   1. cndplottopo (preferred, if available)
%   2. topoplot (standard EEGLAB function)
%   3. bar plot (ultimate fallback)

% Clear the axes
axes(ax);
cla;

% Try cndplottopo first (from CNSP toolbox)
try
    cndplottopo(topoData, chanlocs, ...
        'title', sprintf('Topography @ %.1f ms', timepoint), ...
        'hover', 'off', ...
        'labels', 'off', ...
        'tempMode', 'on');
    colorbar;
    return;
catch ME1
    % cndplottopo failed, try topoplot
    fprintf('Note: cndplottopo failed (%s), using topoplot\n', ME1.message);
end

% Try topoplot (standard EEGLAB function)
try
    topoplot(topoData, chanlocs, ...
        'style', 'both', ...
        'electrodes', 'on', ...
        'plotrad', 0.55, ...
        'headrad', 0.5, ...
        'maplimits', 'maxmin');
    title(sprintf('Topography @ %.1f ms', timepoint), 'FontSize', 12);
    colorbar;
    return;
catch ME2
    % topoplot failed, use bar plot as fallback
    fprintf('Note: topoplot failed (%s), using bar plot\n', ME2.message);
end

% Ultimate fallback: bar plot
bar(topoData);
xlabel('Channel Index', 'FontSize', 11);
ylabel('TRF Weight', 'FontSize', 11);
title(sprintf('Channel Values @ %.1f ms (bar plot)', timepoint), 'FontSize', 12);
grid on;

end

