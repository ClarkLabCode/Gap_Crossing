% Get current figure and axes from tiledlayout
srcFig = gcf;
srcLayout = findobj(srcFig, 'Type', 'tiledlayout');
srcAxes = flipud(findobj(srcLayout, 'Type', 'axes'));

if numel(srcAxes) < 2
    error('Not enough tiles. Need at least two axes.');
end

% Create new figure with 1 row, 3 columns
newFig = figure;
newLayout = tiledlayout(newFig, 1, 3);  % 1 row, 3 columns

%% Tile 1: spans 2 columns (wide)
destAx1 = nexttile(newLayout, [1 2]);  % span 2 columns
srcAx1 = srcAxes(1);
ax1Children = allchild(srcAx1);
plotted = false;

for obj = flipud(ax1Children')
    if isa(obj, 'matlab.graphics.chart.primitive.Line')
        x = obj.XData;
        y = obj.YData;

        factor = 4;
        x_ds = x(1:factor:end);
        y_ds = y(1:factor:end);

        plot(destAx1, x_ds, y_ds, 'DisplayName', obj.DisplayName);
        hold(destAx1, 'on');
        plotted = true;
    end
end

copy_axis_format(srcAx1, destAx1);

if ~plotted
    title(destAx1, 'Tile 1: No line data to downsample');
else
    title(destAx1, 'Tile 1: Downsampled by 4');
end

%% Tile 2: occupies last column (normal width)
destAx2 = nexttile(newLayout, 3);
srcAx2 = srcAxes(4);
ax2Children = allchild(srcAx2);

for obj = flipud(ax2Children')
    try
        copyobj(obj, destAx2);
    catch
        warning("Skipped unsupported object of class: %s", class(obj));
    end
end

copy_axis_format(srcAx2, destAx2);
title(destAx2, 'Tile 2: Original Copy');

%% Helper function to copy formatting
function copy_axis_format(srcAx, destAx)
    try, destAx.XLim = srcAx.XLim; end
    try, destAx.YLim = srcAx.YLim; end
    try, destAx.XScale = srcAx.XScale; end
    try, destAx.YScale = srcAx.YScale; end
    try, destAx.XTick = srcAx.XTick; end
    try, destAx.YTick = srcAx.YTick; end
    try, destAx.XTickLabel = srcAx.XTickLabel; end
    try, destAx.YTickLabel = srcAx.YTickLabel; end
    try, destAx.Box = srcAx.Box; end
    try, destAx.FontSize = srcAx.FontSize; end

    try, destAx.XLabel.String = srcAx.XLabel.String; end
    try, destAx.YLabel.String = srcAx.YLabel.String; end
    try, destAx.XLabel.FontSize = srcAx.XLabel.FontSize; end
    try, destAx.YLabel.FontSize = srcAx.YLabel.FontSize; end
end
