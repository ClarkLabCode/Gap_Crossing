% Plot the crossing rate by taking all crossings and multiplying by the
% approximate fraction of non-glass crossings at each gap size since NN  
% does a poor job distinguishing between them. This rate was gathered by
% hand labeling more than 300 crossing events at random from the horizontal
% gap crossing experiments.
shadedErrorBar(WS.crossStats.gapSizes,mean(WS.crossStats.AllAllVecAllCrossOverAllGapEventRate,1).*fracCrossing',WS.crossStats.AllAllVecAllCrossOverAllGapEventStd,'lineprops', {'Color','k','Marker','none'})
ylim([0,1])
xticks([1:0.5:2.5])
xlim([0.8,2.7])
xlabel('Gap Width (mm)');
ylabel('Crossings / All Gap Events');