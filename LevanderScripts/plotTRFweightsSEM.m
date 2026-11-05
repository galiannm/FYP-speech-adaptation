function plotTRFweightsSEM(models, chanlocs, semFlag, channels, plotNormFlag, fontSize, thicknessScale, densityType)
%PLOTTRFWEIGHTSSEM Plot TRF weights with multiple visualization options
%
% Inputs:
%   models      - Struct array of TRF models, one per subject
%   chanlocs    - EEGLAB-style channel locations struct
%   semFlag     - 0: mean only, 1: mean+SEM shading, 2: density line
%   channels    - Cell array of channel labels (empty = all)
%   plotNormFlag - Normalization flag (default: 1)
%   fontSize    - Font size (default: 20)
%   thicknessScale - Max line thickness (default: 3)
%   densityType - 'confidence' or 'spread' (default: 'spread')

    % Set defaults
    if nargin < 8 || isempty(densityType)
        densityType = 'spread';
    end
    if nargin < 7 || isempty(thicknessScale)
        thicknessScale = 3;
    end
    if nargin < 6 || isempty(fontSize)
        fontSize = 20;
    end
    if nargin < 5 || isempty(plotNormFlag)
        plotNormFlag = 1;
    end
    if nargin < 4
        channels = [];
    end
    if nargin < 3
        semFlag = 0;
    end

    nSubs = numel(models);

    % Compute model average
    avgModel = mTRFmodelAvg(models, plotNormFlag);
    t = avgModel.t(:)';
    Wmean = squeeze(avgModel.w);
    
    if size(Wmean,2) ~= numel(t)
        Wmean = Wmean';
    end
    [nCh, nT] = size(Wmean);

    % Select channels by label
    allLabels = {chanlocs.labels};
    if isempty(channels)
        chanIdx = 1:nCh;
        chanLabels = allLabels;
        showLegend = false;
    else
        chanIdx = find(ismember(allLabels, channels));
        if isempty(chanIdx)
            error('None of the specified channel labels were found.');
        end
        chanLabels = allLabels(chanIdx);
        showLegend = true;
    end

    % Compute SEM if needed
    if semFlag == 1 || semFlag == 2
        Wall = nan(nSubs, nCh, nT);
        for s = 1:nSubs
            tmp = mTRFmodelAvg(models(s), plotNormFlag);
            Ws = squeeze(tmp.w);
            if size(Ws,2) ~= numel(tmp.t)
                Ws = Ws';
            end
            minLen = min(nT, size(Ws,2));
            Wall(s,:,1:minLen) = Ws(:,1:minLen);
        end

        validLen = find(all(~isnan(Wall(:,1,:)),1),1,'last');
        t = t(1:validLen);
        Wmean = Wmean(:,1:validLen);
        Wall = Wall(:,:,1:validLen);

        % Compute SEM
        Wsem = squeeze(std(Wall,0,1)) ./ sqrt(nSubs);
    end

    % DO NOT create figure - plot into current axes
    hold on;
    
    % Plot zero reference lines (bolder)
    xline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    cols = lines(numel(chanIdx));
    legendItems = {};
    legendHandles = [];
    
    % Create dummy handles for legend items
    if semFlag == 1
        hSEM = fill(NaN, NaN, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
        legendItems{end+1} = '±SEM';
        legendHandles(end+1) = hSEM;
        
        hTRF = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
        legendItems{end+1} = 'TRF Weight';
        legendHandles(end+1) = hTRF;
    elseif semFlag == 2
        hDensity = plot(NaN, NaN, 'k-', 'LineWidth', 2);
        legendItems{end+1} = 'TRF Weight (density)';
        legendHandles(end+1) = hDensity;
    else
        hMean = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
        legendItems{end+1} = 'TRF Weight';
        legendHandles(end+1) = hMean;
    end

    % Plot actual data for each channel
    for ii = 1:numel(chanIdx)
        ch = chanIdx(ii);
        label = chanLabels{ii};
        yMean = Wmean(ch,:);
        
        % SEM shading
        if semFlag == 1
            ySem = Wsem(ch,:);
            upper = yMean + ySem;
            lower = yMean - ySem;
            
            fill([t fliplr(t)], [upper fliplr(lower)], ...
                 cols(ii,:), 'EdgeColor', 'none', 'FaceAlpha', 0.25, ...
                 'HandleVisibility', 'off');
        end
        
        % Density line
        if semFlag == 2
            ySem = Wsem(ch,:);
            
            if strcmpi(densityType, 'confidence')
                semInv = 1./(ySem + eps);
                lineWidths = 1 + thicknessScale * ...
                            (semInv - min(semInv)) ./ (max(semInv) - min(semInv));
            else
                lineWidths = 1 + thicknessScale * ...
                            (ySem - min(ySem)) ./ (max(ySem) - min(ySem));
            end
            
            % Plot line with varying thickness
            for seg = 1:length(t)-1
                plot(t(seg:seg+1), yMean(seg:seg+1), ...
                     'Color', cols(ii,:), ...
                     'LineWidth', lineWidths(seg), ...
                     'HandleVisibility', 'off');
            end
            hChan = plot(NaN, NaN, '-', 'Color', cols(ii,:), 'LineWidth', 2);
        else
            hChan = plot(t, yMean, '-', 'Color', cols(ii,:), 'LineWidth', 1.5);
        end
        
        % Add channel to legend
        legendItems{end+1} = label;
        legendHandles(end+1) = hChan;
    end

    % Format plot
    xlabel('Time-latency (ms)', 'FontSize', fontSize)
    ylabel('TRF weights (a.u.)', 'FontSize', fontSize)
    
    % Title
    if showLegend
        if isscalar(chanIdx)
            title(['TRF: ', chanLabels{1}], 'FontSize', fontSize + 2, 'FontWeight', 'bold');
        else
            if numel(chanLabels) <= 4
                chanStr = strjoin(chanLabels, ', ');
            else
                chanStr = [strjoin(chanLabels(1:3), ', '), ', ...'];
            end
            title(['TRF Weights: ', chanStr], 'FontSize', fontSize + 2, 'FontWeight', 'bold');
        end
    else
        title('TRF Weights (All Channels)', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
    end
    
    grid on
    set(gca, 'FontSize', fontSize)
    set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '-');
    
    % Add legend
    if showLegend || semFlag > 0
        hLeg = legend(legendHandles, legendItems, 'Location', 'best');
        set(hLeg, 'FontSize', fontSize - 2, 'Box', 'on', 'EdgeColor', [0.3 0.3 0.3]);
    end
    
    hold off;
end