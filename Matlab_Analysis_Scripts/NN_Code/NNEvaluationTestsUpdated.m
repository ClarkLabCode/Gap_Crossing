% Must be used with the NN Training workspace loaded in (from 10/3/22)

% Compute everything needed for the subfigures corresponding to the neural
% net that classifies crossings vs glass circumventions

% Grab the number and fraction of cases labeled as ambiguous
numTrueAmbig = sum(imdsTest_AmbUnamb_0_1.Labels == 'ambig');
fracTrueAmbig = numTrueAmbig/length(imdsTest_AmbUnamb_0_1.Labels);

% Have the neural net output its predictions on the test data
PredictionsOfNetWithAmbig = predict(net_OnOffAmbig, imdsTest_OnOff);

% Set the different ambiguity thresholds to sweep through
thresh = [0.01, 0.025, 0.05, 0.1:0.1:0.3, 1/3, 0.4:0.1:0.6, 2/3, 0.7:0.1:0.9, 0.95, 0.975, 0.99];

% Initialize vectors to hold relevant metrics for each ambiguity threshold
fracOfDataUnambig = zeros(length(thresh),1);
AUC =  zeros(length(thresh),1);

% Grab the number of cases that were labeled off glass, on glass, and ambig
numOffGlass = sum(imdsTest_OnOff.Labels == 'off_glass_crossing');
numOnGlass = sum(imdsTest_OnOff.Labels == 'on_glass_crossing');
numAmbig = sum(imdsTest_AmbUnamb_0_1.Labels == 'ambig');

% Initialize vectors that hold which of the test imds has which labels
imgNumOffGlass = zeros(numOffGlass,1);
imgNumOnGlass = zeros(numOnGlass,1);
imgNumAmbig = zeros(numAmbig,1);

% Go through and fill the vectors initialized above for off glass, on
% glass, and ambig
% Note that all the imds file names have a 5 digit number associated with
% them that is independent of their labels, so this is what's being grabbed
% by the (end-8:end-4) bit in the loops below
for i = 1:numOffGlass
    imgNumOffGlass(i) = str2double(imdsTest_OnOff.Files{i}(end-8:end-4));
end
for i = 1:numOnGlass
    imgNumOnGlass(i) = str2double(imdsTest_OnOff.Files{i+numOffGlass}(end-8:end-4));
end
for i = 1:numAmbig
    imgNumAmbig(i) = str2double(imdsTest_AmbUnamb_0_1.Files{i}(end-8:end-4));
end

% Now compute how many of the off glass and on glass cases were also
% labeled as ambiguous cases
numOffGlassAmbigToo = length(intersect(imgNumAmbig,imgNumOffGlass));
numOnGlassAmbigToo = length(intersect(imgNumAmbig,imgNumOnGlass));

%% Generate ROC for various threshold values of ambiguity

figure

% Go through each different ambiguity threshold and plot the ROC for the
% cases that were not labeled as ambiguous by the neural net
hold on
for i = 1:length(thresh)
    % Determine which cases are labeled as ambiguous by the neural net
    AmbigCases = PredictionsOfNetWithAmbig(:,1) > thresh(i);
    % Generate the ROC and AUC among the non-ambiguous cases
    [x,y,~,auc] = perfcurve(imdsLabelsMat(~AmbigCases,1)', PredictionsOfNetWithAmbig(~AmbigCases,2)',1);
    % Plot the ROC and use a red to blue color gradient for each curve
    plot(x,y, 'LineWidth', 2, 'color', [1-(i-1)/(length(thresh)-1) 0 (i-1)/(length(thresh)-1)])
    % Compute the fraction of data that the neural net labeled as unambig
    fracOfDataUnambig(i) = sum(~AmbigCases)/length(AmbigCases);
    % Fill the AUC vector with the AUC from each individual threshold
    AUC(i) = auc;
end

% Below is a bunch of coding gymnastics to allow the plot with all the ROCs
% to have a legend that shows the ambiguity threshold, the fraction of data
% that is below that ambiguity threshold (in parentheses as a %), and the
% ROC of the neural net that isn't trained with ambiguity
openParenth = '      (';
closeParenth = '%)';
numericLabels = ...
    [num2str(thresh','%.3f'), ...
     convertCharsToStrings(cellstr(repmat(openParenth,[length(thresh),1]))), ...
     erase(convertCharsToStrings(cellstr(char(num2str(round(100*fracOfDataUnambig))))),' '),...
     convertCharsToStrings(cellstr(repmat(closeParenth,[length(thresh),1])))];
numericLabels = append(numericLabels(:,1),numericLabels(:,2),numericLabels(:,3),numericLabels(:,4));
legendLabels = char(char(convertStringsToChars(numericLabels)),'Not Trained for Ambiguity');

% Now grab the stuff we need for plotting the ROC of the neural net not 
% trained with ambiguity
imdsLabelsMat = zeros(200,2);
imdsLabelsMat(:,1) = (imdsTest_OnOff.Labels == 'off_glass_crossing');
imdsLabelsMat(:,2) = (imdsTest_OnOff.Labels == 'on_glass_crossing');

% Plot the ROC for the neural net not trained with ambiguity in black
[x,y,~,AUC_No_Ambig] = perfcurve(imdsLabelsMat(:,1)', PredictionsOfNetOnOff(:,1)',1);
plot(x,y,'k','LineWidth',2);

% Add the axes labels, title, and legend for the ROCs
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC for Varying Ambiguity Thresholds');
legend(cellstr(legendLabels));
hold off

%% Calculate #Off/(#On+#Off) as a function of ambiguity thresh

% Initialize a vector holding P(crossing) as a function of ambig thresh
FracOffOverOn = zeros(length(thresh),1);

% Go through each ambig thresh and compute P(crossing) for that thresh by
% counting the number of unambiguous cases that were labeled as off glass
for i = 1:length(thresh)
    % Determine which cases are labeled as ambiguous by the neural net
    AmbigCases = PredictionsOfNetWithAmbig(:,1) > thresh(i);
    % Among non-ambiguous cases, see which ones were labeled as off glass
    % by the neural net
    OffGlassPredAmongUnambig = PredictionsOfNetWithAmbig(~AmbigCases,2)>PredictionsOfNetWithAmbig(~AmbigCases,3);
    % Compute the fraction of off glass crossings (i.e., P(crossing))
    % Note that the denominator is all the cases, not just the unambiguous
    % cases, because we aren't disposing the ambiguous cases but rather
    % labeling them as on glass (i.e., glass circumventions)
    FracOffOverOn(i) = sum(OffGlassPredAmongUnambig)/(length(AmbigCases));
end

% Now we go ahead and plot everything that we just collected above
figure

% Plot P(crossing) as a function of ambiguity threshold, matching the color
% gradient in the previous figure from red to blue
hold on
for i = 1:(length(thresh)-1)
    plot(thresh(i:i+1), FracOffOverOn(i:i+1),'color',[1-(i-1)/(length(thresh)-1) 0 (i-1)/(length(thresh)-1)], 'LineWidth', 2)
end

% Now add a black dashed line showing the ratio for the true labels
yline(sum(imdsLabelsMat(:,1))/length(imdsLabelsMat),'k--', 'LineWidth', 2);

% Finally, add a gray error bar around the true labels line which shows the
% maximal uncertainty in the true labels based on what were hand labeled as
% ambiguous cases
x_vec = 0:0.1:1;
shadedErrorBar(x_vec,...
               repmat(sum(imdsLabelsMat(:,1))/length(imdsLabelsMat),[length(x_vec),1]),... % true labels
              [repmat(numOnGlassAmbigToo/length(imdsLabelsMat),[length(x_vec),1]),... % upper bound
               repmat(numOffGlassAmbigToo/length(imdsLabelsMat),[length(x_vec),1])],... % lower bound
               'lineProps',{'LineStyle','none','Color','k'})
hold off

% Now just add the axes labels
xlabel('Ambiguity Threshold');
ylabel('P(Crossing)');
xlim([0 1]);
ylim([0 1]);

%% Plot AUC as a function of ambiguity threshold

figure

% Go through and plot the AUC as a function of ambiguity thresh with the
% same color scheme as before from red to blue
hold on
for i = 1:(length(thresh)-1)
    plot(thresh(i:i+1), AUC(i:i+1),'color',[1-(i-1)/(length(thresh)-1) 0 (i-1)/(length(thresh)-1)], 'LineWidth', 2)
end

% Plot the AUC for the neural net that was not trained with ambiguity
yline(AUC_No_Ambig,'k-','LineWidth', 2)
hold off

% Add the axes labels
xlabel('Ambiguity Threshold');
ylabel('AUC within Non-Ambiguous Cases');
xlim([0 1]);
ylim([0 1]);


%% Generate ROC for each different gap size

% Initialize a vector that holds the AUC for each gap size
AUC_gap = zeros(NumGaps,1);

figure

% Go through each gap size and generate the ROC
hold on
for gapCounter = 1:NumGaps
    % First only look at the non-ambiguous cases
    AmbigCases = PredictionsOfNetWithAmbig(:,1) > 0.5;
    % Find the ones corresponding to the right gap number
    YesToGapID = (gapIDSampleVec == gapCounter);
    % Only keep non-ambiguous cases with the right gap number
    EventOfInterest = (~AmbigCases & YesToGapID);
    % Generate ROC
    [x,y,~,auc] = perfcurve(imdsLabelsMat(EventOfInterest,1)', PredictionsOfNetWithAmbig(EventOfInterest,2)',1);
    % Plot it
    plot(x,y, 'LineWidth', 2)
    % Save the AUC
    AUC_gap(gapCounter) = auc;
end

% Add the axes labels, title, and legend
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC for Each Gap Size');
legend({'1 mm Gap'; '1.5 mm Gap'; '2 mm Gap'; '2.5 mm Gap'});
hold off