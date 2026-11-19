function [violinFig, contributorCounts] = plotTrialViolins(n1Data,p2Data,corrData,seriesLabels)
%PLOTTRIALVIOLINS Descriptive distributions across series at each trial position.
% Rows are dataset/feature series; columns are presentation positions.
% Colours identify series. Black bars mark medians. Non-finite values are omitted.
% A Gaussian density is drawn only for at least two distinct finite values.

assert(isequal(size(n1Data),size(p2Data),size(corrData)), 'N1, P2 and correlation arrays must have matching dimensions.');
assert(size(n1Data,1) == numel(seriesLabels), 'Provide one label per measurement row.');

allData = {n1Data,p2Data,corrData};
panelTitles = {'N1 magnitude','P2 magnitude','Prediction correlation'};
yLabels = {'Normalised N1 GFP maximum', 'Normalised P2 GFP maximum','95th-percentile-scaled prediction correlation'};
nSeries = size(n1Data,1);
hasAnyData = false(nSeries,1);
lastTrial = 0;
for panel = 1:3
    hasAnyData = hasAnyData | any(isfinite(allData{panel}),2);
    cols = find(any(isfinite(allData{panel}),1));
    if ~isempty(cols)
        lastTrial = max(lastTrial,cols(end));
    end
end
if lastTrial == 0
    warning('plotTrialViolins:NoData','No finite measurements available for violin plots.');
    violinFig = gobjects(0);
    contributorCounts = zeros(3,0);
    return
end

colours = lines(nSeries);
if nSeries > 7
    colours = turbo(nSeries);
end
% Fixed horizontal offsets keep each series identifiable without changing RNG.
if nSeries == 1
    offsets = 0;
else
    offsets = linspace(-0.15,0.15,nSeries);
end

violinFig = figure('Name','Trial-by-trial envelope violin summaries', ...
    'Units','normalized','Position',[0.05,0.08,0.90,0.80]);
violinLayout = tiledlayout(violinFig,3,1, ...
    'TileSpacing','compact','Padding','compact');
title(violinLayout, ...
    'Trial-by-trial envelope analysis: dots = series; black bars = medians');
contributorCounts = zeros(3,lastTrial);

for panel = 1:3
    ax = nexttile(violinLayout);
    hold(ax,'on');
    data = allData{panel};
    for trial = 1:lastTrial
        valid = isfinite(data(:,trial));
        values = data(valid,trial);
        n = numel(values);
        contributorCounts(panel,trial) = n;
        if n == 0
            continue
        end

        if n >= 2 && std(values) > 0
            bandwidth = 1.06*std(values)*n^(-1/5);
            densityY = linspace(min(values)-3*bandwidth, ...
                max(values)+3*bandwidth,100)';
            z = (densityY-values')/bandwidth;
            density = mean(exp(-0.5*z.^2),2)/(bandwidth*sqrt(2*pi));
            halfWidth = 0.34*density/max(density);
            fill(ax,[trial-halfWidth;flipud(trial+halfWidth)], ...
                [densityY;flipud(densityY)],[0.72,0.78,0.90], ...
                'FaceAlpha',0.35,'EdgeColor',[0.45,0.50,0.65], ...
                'HandleVisibility','off');
        end

        scatter(ax,trial+offsets(valid),values,18,colours(valid,:), ...
            'filled','HandleVisibility','off');
        plot(ax,trial+[-0.18,0.18],[median(values),median(values)], ...
            'k-','LineWidth',1.2,'HandleVisibility','off');
    end

    if panel < 3
        yline(ax,0,'k:','HandleVisibility','off');
    end
    xlim(ax,[0.5,lastTrial+0.5]);
    xticks(ax,1:lastTrial);
    tickLabels = arrayfun(@(j) sprintf('%d (n=%d)', ...
        j,contributorCounts(panel,j)),1:lastTrial,'UniformOutput',false);
    xticklabels(ax,tickLabels);
    xtickangle(ax,45);
    xlabel(ax,'Trial presentation position (n = contributing series)');
    ylabel(ax,yLabels{panel});
    title(ax,panelTitles{panel});
    grid(ax,'on');

    if panel == 1
        legendHandles = gobjects(sum(hasAnyData),1);
        rows = find(hasAnyData);
        for k = 1:numel(rows)
            legendHandles(k) = plot(ax,NaN,NaN,'o', ...
                'MarkerFaceColor',colours(rows(k),:), ...
                'MarkerEdgeColor',colours(rows(k),:),'LineStyle','none');
        end
        lgd = legend(ax,legendHandles,seriesLabels(hasAnyData), ...
            'Interpreter','none','FontSize',8);
        lgd.Layout.Tile = 'east';
    end
end
end
