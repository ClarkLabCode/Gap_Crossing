% === BEGIN SCRIPT ===
% 1) Grab all axes (including those inside a tiledlayout)
fig       = gcf;
axHandles = findall(fig, 'Type', 'Axes');
% 2) Pull out every X/Y from each child object *except* red lines
allX = [];
allY = [];
for k = 1:numel(axHandles)
ch = axHandles(k).Children;
for c = 1:numel(ch)
obj = ch(c);
% skip any pure‐line object drawn in red (per‐tile fit lines)
if isprop(obj, 'Color') && isequal(obj.Color, [1 0 0])
continue;
end
% Only take things that actually have XData/YData
if isprop(obj, 'XData') && isprop(obj, 'YData')
x = obj.XData(:);
y = obj.YData(:);
allX = [allX; x];
allY = [allY; y];
end
end
end
% 3) Scatter them all together
figure;
scatter(allX(allX<1.5), allY(allX<1.5), 'blue');
hold on;
scatter(allX(allX>=1.5), allY(allX>=1.5), 'red');
% 4) Fit linear model (requires Statistics Toolbox)
lm = fitlm(allX, allY);
% 5) Extract slope, intercept, t‐stat and df
intercept = lm.Coefficients.Estimate(1);
slope     = lm.Coefficients.Estimate(2);
tStat     = lm.Coefficients.tStat(2);
df        = lm.DFE;
% 6) Compute one‐sided p‐value for H1: slope < 0
pOneSide = tcdf(tStat, df);
% 7) Plot overall best‐fit line in red
xfit = linspace(min(allX), max(allX), 100);
yfit = slope*xfit + intercept;
plot(xfit, yfit, 'g-', 'LineWidth', 2);
% 8) Annotate slope and p‐value on the plot
slopeStr = sprintf('Slope = %.3f', slope);
pStr     = sprintf('p = %.2e', pOneSide);  % scientific notation, 3 sig figs
% Choose a location: 5% in from left, 10% down from top
xlims = xlim;
ylims = ylim;
ylim([-1, 10])
yticks(0:5:10)
yticklabels(0:5:10)
xpos  = xlims(1) + 0.05*(xlims(2)-xlims(1));
ypos  = ylims(2) - 0.10*(ylims(2)-ylims(1));
text(xpos, ypos, {slopeStr, pStr}, ...
'FontSize', 12, ...
'BackgroundColor', 'white', ...
'EdgeColor', 'black');
% 9) Tidy up
xlabel('X');
ylabel('Y');
title('Combined Scatter with Overall Best-Fit');
legend({'Data','Overall Fit'}, 'Location','best');
hold off;
% 10) Also print to console
fprintf('Estimated slope      = %.3f\n', slope);
fprintf('One-sided p-value   = %.2e\n', pOneSide);
% === END SCRIPT ===