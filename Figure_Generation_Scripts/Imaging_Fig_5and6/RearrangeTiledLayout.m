% Get current figure and layout
fig = gcf;
oldLayout = findall(fig, 'Type', 'tiledlayout');
oldAxes = findall(oldLayout, 'Type', 'axes');

% Get axes in the correct visual order (reverse of findall)
oldAxes = flipud(oldAxes);

% Create new layout (change these values as desired)
newRows = 3;
newCols = 6;
newLayout = tiledlayout(fig, newRows, newCols);

% Copy plots into new layout
for k = 1:length(oldAxes)
    oldAx = oldAxes(k);
    newAx = nexttile(newLayout);

    % Copy all children (lines, images, etc.)
    for child = flipud(allchild(oldAx))'
        copyobj(child, newAx);
    end

    % Copy labels, limits, titles, etc.
    newAx.XLim = oldAx.XLim;
    newAx.YLim = oldAx.YLim;
    newAx.XLabel.String = oldAx.XLabel.String;
    newAx.YLabel.String = oldAx.YLabel.String;
    newAx.Title.String = oldAx.Title.String;

    % Optional: copy axis ticks, scales, etc.
    newAx.XTick = oldAx.XTick;
    newAx.YTick = oldAx.YTick;
    newAx.XScale = oldAx.XScale;
    newAx.YScale = oldAx.YScale;
end

% Optionally delete old layout
delete(oldLayout);
