colorCross = [0,225/255,0];
colorGlass = [0,1,1];
colorCirc = [0,0,1];
colorRet = [1,69/255,0];

colororder([colorCross;colorGlass;colorCirc;colorRet]);

% numPlots = 10;
numPlots = 1;

% TL = tiledlayout(2,5);
allWSPlotted = cell(numPlots,1);
allWSPlotted{1} = WS_WT;
% allWSPlotted{6} = WS_Dark;
% allWSPlotted{2} = WS_One_Eye;
% allWSPlotted{7} = WS_Two_Eye;
% allWSPlotted{8} = WS_HF;
% allWSPlotted{9} = WS_T4T5_shR;
% allWSPlotted{10} = WS_T4T5_dTRPA1;
% allWSPlotted{3} = WS_T4T5_ID1;
% allWSPlotted{4} = WS_sp_E_shR;
% allWSPlotted{5} = WS_sp_E_dTRPA1;

plotTitles = cell(numPlots,1);
plotTitles{1} = 'IsoD1 (control)';
% plotTitles{6} = 'IsoD1 (dark)';
% plotTitles{2} = 'IsoD1 (one eye painted)';
% plotTitles{7} = 'IsoD1 (both eyes painted)';
% plotTitles{8} = 'IsoD1 (head-fixed)';
% plotTitles{9} = 'T4T5 > shibire';
% plotTitles{10} = 'T4T5 > dTRPA1';
% plotTitles{3} = 'T4T5 > IsoD1';
% plotTitles{4} = 'sp empty DBD > shibire';
% plotTitles{5} = 'sp empty DBD > dTRPA1';

for j = 1:numPlots
% ax = nexttile;
tempCrossStats = allWSPlotted{j}.crossStats;

activeFlies = tempCrossStats.AllUpVecAllGapEventsLog;

AllGapEventFreqs = zeros(4,4);
AllGapEventFreqsSEM = zeros(4,4);
gapSizes = 1:0.5:2.5;

AllGapEventFreqs(:,1) = mean(tempCrossStats.AllUpVecProperCross(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1);
AllGapEventFreqs(:,2) = mean(tempCrossStats.AllUpVecGlassCross(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1);
AllGapEventFreqs(:,3) = mean(tempCrossStats.AllUpVecCircumvent(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1);
AllGapEventFreqs(:,4) = mean(tempCrossStats.AllUpVecRetreat(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1);

AllGapEventFreqsSEM(:,1) = std(tempCrossStats.AllUpVecProperCross(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1)/sqrt(size(tempCrossStats.AllUpVecProperCross(activeFlies,:),1));
AllGapEventFreqsSEM(:,2) = std(tempCrossStats.AllUpVecGlassCross(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1)/sqrt(size(tempCrossStats.AllUpVecProperCross(activeFlies,:),1));
AllGapEventFreqsSEM(:,3) = std(tempCrossStats.AllUpVecCircumvent(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1)/sqrt(size(tempCrossStats.AllUpVecProperCross(activeFlies,:),1));
AllGapEventFreqsSEM(:,4) = std(tempCrossStats.AllUpVecRetreat(activeFlies,:)./tempCrossStats.AllUpVecAllGapEvents(activeFlies,:),1)/sqrt(size(tempCrossStats.AllUpVecProperCross(activeFlies,:),1));

bar(gapSizes,AllGapEventFreqs,'stacked')
hold on
yVal = zeros(4,1);
for i = 1:3
    yVal = yVal + AllGapEventFreqs(:,i);
    er = errorbar(gapSizes,yVal,AllGapEventFreqsSEM(:,i));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
end
hold off

% xlabel('Gap Size (mm)')
% ylabel('Frequency');
title(plotTitles{j})
xlim([0.6,2.9])
ylim([0,1.1])
xticks(gapSizes)
yticks(0:0.5:1)
xlabel('Gap Width (mm)');
ylabel('Frequency');
end

lg = legend({'Cross','Glass Circ','Gap Circ','Retreat'});

                lg.FontSize = 10;
                        lg.Layout.Tile = 'North'; % <-- Legend placement with tiled layout
                        lg.Orientation = 'horizontal';

                % Add a global x and y axis title since all experiments share those axes
% TL.XLabel.String = 'Gap Width (mm)';
% TL.YLabel.String = 'Frequency';