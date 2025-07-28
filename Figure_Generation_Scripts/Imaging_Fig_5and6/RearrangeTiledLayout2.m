% Get current figure and tiledlayout
srcFig = gcf;
srcLayout = findobj(srcFig, 'Type', 'tiledlayout');
srcAxes = findobj(srcLayout, 'Type', 'axes');

% Filter valid axes only and reverse to preserve order
srcAxes = flipud(srcAxes(isgraphics(srcAxes, 'axes')));

% Define new layout shape
newRows = 3;                  % <--- Change this
newCols = 6;    % or set to e.g., 2 for 2x2, etc.

% Create new figure with new tiledlayout
destFig = figure;
destLayout = tiledlayout(destFig, newRows, newCols);

% Copy each axes into new tiledlayout
for k = 1:numel(srcAxes)
    srcAx = srcAxes(k);
    destAx = nexttile(destLayout);

    % Copy graphical objects (lines, bars, images, etc.)
    children = allchild(srcAx);
    copyobj(children, destAx);

    % Copy axis limits and labels
    destAx.XLim = srcAx.XLim;
    destAx.YLim = srcAx.YLim;
    destAx.XScale = srcAx.XScale;
    destAx.YScale = srcAx.YScale;

    destAx.XLabel.String = srcAx.XLabel.String;
    destAx.YLabel.String = srcAx.YLabel.String;
    destAx.Title.String  = srcAx.Title.String;

    % Optional: Match ticks, grid, etc.
    destAx.XTick = srcAx.XTick;
    destAx.YTick = srcAx.YTick;
    destAx.XTickLabel = srcAx.XTickLabel;
    destAx.YTickLabel = srcAx.YTickLabel;
    destAx.Box = srcAx.Box;
    destAx.FontSize = srcAx.FontSize;
end
