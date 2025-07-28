x_thresh = 1.5;

fig = gcf;
layout = findobj(fig, 'Type', 'tiledlayout');
axesHandles = flipud(findobj(layout, 'Type', 'axes'));

for axIdx = 1:numel(axesHandles)
    ax = axesHandles(axIdx);
    hold(ax, 'on');

    % Handle scatter plots (color by x_thresh)
    scatters = findobj(ax, 'Type', 'scatter');
    for sc = scatters'
        x = sc.XData;
        y = sc.YData;

        idx_below = x < x_thresh;
        idx_above = x >= x_thresh;

        markerSize = sc.SizeData;
        if numel(markerSize) == 1
            markerSize = repmat(markerSize, size(x));
        end

        % Plot blue points below threshold
        if any(idx_below)
            scatter(ax, x(idx_below), y(idx_below), markerSize(idx_below), 'b');
        end

        % Plot red points at/above threshold
        if any(idx_above)
            scatter(ax, x(idx_above), y(idx_above), markerSize(idx_above), 'r');
        end

        delete(sc);  % Remove original scatter
    end

    % Handle line plots (make entire line green)
    lines = findobj(ax, 'Type', 'line');
    for lineObj = lines'
        x = lineObj.XData;
        y = lineObj.YData;
        lw = lineObj.LineWidth;

        plot(ax, x, y, 'Color', [0 0.7 0], 'LineWidth', lw);  % green line

        delete(lineObj);  % Remove original line
    end
end
