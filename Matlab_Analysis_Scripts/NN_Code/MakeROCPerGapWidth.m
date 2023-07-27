AUC_gap = zeros(NumGaps,1);

figure
hold on
for i = 1:NumGaps
    AmbigCases = PredictionsOfNetWithAmbig(:,1) > 0.5;
    YesToGapID = (gapIDSampleVec == i);
    EventOfInterest = (~AmbigCases | YesToGapID);
    [x,y,~,auc] = perfcurve(imdsLabelsMat(EventOfInterest,1)', PredictionsOfNetWithAmbig(EventOfInterest,2)',1);
    plot(x,y, 'LineWidth', 2)
    AUC_gap(i) = auc;
end
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve for Each Gap Size');
legend({'1 mm Gap'; '1.5 mm Gap'; '2 mm Gap'; '2.5 mm Gap'});
hold off