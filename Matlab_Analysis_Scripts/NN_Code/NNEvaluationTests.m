%% Generate ROC curves for various threshold values of ambiguity

PredictionsOfNetWithAmbig = predict(net_OnOffAmbig, imdsTest_OnOff);

thresh = [0.01:0.01:0.09, 0.1:0.1:0.9, 0.95, 0.975, 0.99];
fracOfDataUnambig = zeros(length(thresh),1);
AUC =  zeros(length(thresh),1);

figure

hold on
for i = 1:length(thresh)
AmbigCases = PredictionsOfNetWithAmbig(:,1) > thresh(i);
% figure(i)
% plotroc(imdsLabelsMat(~AmbigCases,1)', PredictionsOfNetWithAmbig(~AmbigCases,2)');
[x,y,~,auc] = perfcurve(imdsLabelsMat(~AmbigCases,1)', PredictionsOfNetWithAmbig(~AmbigCases,2)',1);
plot(x,y, 'LineWidth', 2)
% title(['ROC Curve for 3 Class Classifier with Ambiguity Threshold ', num2str(thresh(i))]);
fracOfDataUnambig(i) = sum(~AmbigCases)/length(AmbigCases);
AUC(i) = auc;
end
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve for Varying Ambiguity Thresholds');
legend(num2str(thresh'));
hold off

%% Calculate #Off/(#On+#Off) among unambiguous cases as a function of ambiguity thresh

PredictionsOfNetWithAmbig = predict(net_OnOffAmbig, imdsTest_OnOff);

thresh = [0.01:0.01:0.09, 0.1:0.1:0.9, 0.95, 0.975, 0.99];
FracOffOverOn = zeros(length(thresh),1);

for i = 1:length(thresh)
AmbigCases = PredictionsOfNetWithAmbig(:,1) > thresh(i);
OffGlassPredAmongUnambig = PredictionsOfNetWithAmbig(~AmbigCases,2)>PredictionsOfNetWithAmbig(~AmbigCases,3);
FracOffOverOn(i) = sum(OffGlassPredAmongUnambig)/(length(AmbigCases)-sum(AmbigCases));
end

figure

plot(thresh, FracOffOverOn)
hold on
errorbar(thresh, sum(imdsLabelsMat(:,1))/length(imdsLabelsMat)*ones(length(thresh),1),0.12*ones(length(thresh),1));
hold off
xlabel('Ambiguity Threshold');
ylabel('Off/(On + Off)');
ylim([0 1]);
legend('Unambiguous Data','All Data (+/- Ambig Data)', 'Location','northeast');

%% Generate ROC curve for classifier trained only with On/Off

imdsLabelsMat = zeros(200,2);
imdsLabelsMat(:,1) = (imdsTest_OnOff.Labels == 'off_glass_crossing');
imdsLabelsMat(:,2) = (imdsTest_OnOff.Labels == 'on_glass_crossing');

figure

% plotroc(imdsLabelsMat(:,1)', PredictionsOfNetOnOff(:,1)');
[x,y,~,auc] = perfcurve(imdsLabelsMat(:,1)', PredictionsOfNetOnOff(:,1)',1);
plot(x,y)
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve for Only On/Off Classifier');

%% Generate ROC curve for each different gap size

AUC_gap = zeros(NumGaps,1);

figure
hold on
for i = 1:NumGaps
    AmbigCases = PredictionsOfNetWithAmbig(:,1) > 0.5;
    YesToGapID = (gapIDSampleVec == i);
    EventOfInterest = (~AmbigCases & YesToGapID);
    [x,y,~,auc] = perfcurve(imdsLabelsMat(EventOfInterest,1)', PredictionsOfNetWithAmbig(EventOfInterest,2)',1);
    plot(x,y, 'LineWidth', 2)
    AUC_gap(i) = auc;
end
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve for Each Gap Size');
legend({'1 mm Gap'; '1.5 mm Gap'; '2 mm Gap'; '2.5 mm Gap'});
hold off

%% Find fraction of Off/(On+Off) for each gap size

FracOffOverOnPerGap = zeros(NumGaps,1);

for i = 1:NumGaps
    AmbigCases = PredictionsOfNetWithAmbig(:,1) > 0.5;
    YesToGapID = (gapIDSampleVec == i);
    EventOfInterest = (~AmbigCases & YesToGapID);
    OffGlassPredForGap = PredictionsOfNetWithAmbig(EventOfInterest,2)>PredictionsOfNetWithAmbig(EventOfInterest,3);
    FracOffOverOnPerGap(i) = sum(imdsLabelsMat(EventOfInterest,1))/sum(EventOfInterest);
end