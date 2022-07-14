WS1 = WS_LC4_shts;
WS2 = WS_LC4_IsoD1;
WS3 = WS_empty_shts;

p1_2 = zeros(4,1);
p1_3 = zeros(4,1);

for i = 1:4
p1_2(i) = ranksum(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i), WS2.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i));
p1_3(i) = ranksum(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i), WS3.crossStats.AllAllVecProperCrossOverAllGapEventRate(:,i));
end

errorbar(1:0.5:2.5,mean(WS1.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS1.crossStats.AllAllVecProperCrossOverAllGapEventStd);
hold on
errorbar((1:0.5:2.5),mean(WS2.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS2.crossStats.AllAllVecProperCrossOverAllGapEventStd);
errorbar((1:0.5:2.5),mean(WS3.crossStats.AllAllVecProperCrossOverAllGapEventRate),WS3.crossStats.AllAllVecProperCrossOverAllGapEventStd);
hold off
ylim([0,1]);
xlim([0.8,2.7]);
xticks(1:0.5:2.5);
xlabel('Gap Width (mm)');
ylabel('Proper Cross Events / All Gap Events');
legend('split LC16 > +; shts; shts', 'split LC16 > +; +; +', 'split empty > +; shts; shts');

[p1_2, p1_3]'