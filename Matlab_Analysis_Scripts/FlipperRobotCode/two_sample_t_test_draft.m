WS1 = WS_LC16_shts;
WS2 = WS_LC16_IsoD1;
WS3 = WS_empty_shts;
WS4 = WS_empty_shts_dark;

p1_2 = zeros(4,1);
p1_3 = zeros(4,1);
p1_4 = zeros(4,1);

for i = 1:4
p1_2(i) = ranksum(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i), WS2.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i));
p1_3(i) = ranksum(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i), WS3.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i));
p1_4(i) = ranksum(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i), WS4.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i));
end

errorbar(1:0.5:2.5,mean(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS1.crossStats.AllAllVecProperCrossOverAllGapEventStd,'Color','r','Marker','.','MarkerSize',15);
hold on
errorbar((1:0.5:2.5),mean(WS2.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS2.crossStats.AllAllVecProperCrossOverAllGapEventStd,'Color',[.7 .7 .7],'Marker','.','MarkerSize',15);
errorbar((1:0.5:2.5),mean(WS3.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS3.crossStats.AllAllVecProperCrossOverAllGapEventStd,'Color','k','Marker','.','MarkerSize',15);
errorbar((1:0.5:2.5),mean(WS4.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS4.crossStats.AllAllVecProperCrossOverAllGapEventStd,'Color','k','LineStyle','--','Marker','.','MarkerSize',15);
hold off
ylim([0,1]);
xlim([0.8,2.7]);
xticks(1:0.5:2.5);
xlabel('Gap Width (mm)');
ylabel('Proper Cross Events / All Gap Events');
legend('LC16 silenced', 'Genetic Control 1', 'Genetic Control 2', 'Genetic Control in Dark');

[p1_2, p1_3, p1_4]'