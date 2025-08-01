allWS = {WS_ID1, WS_ID1Dark, WS_EmptyDBDShr, WS_T4T5shR, WS_T4T5Plus, WS_LC15shR, WS_LC15Plus, WS_EmptyShrDark};
numWS = length(allWS);

MeanXVec = NaN(28,numWS);
StdXVec = NaN(28,numWS);

avgMeanXVec = NaN(numWS, 1);
semMeanXVec = NaN(numWS, 1);

XYTStruct = cell(numWS,2);

for WScounter = 1:numWS
    WS = allWS{WScounter};

    numFlies = length(WS.FlipBinnedFlyStruct);
    numFlips = 0;
    for flyNum = 1:numFlies
        numFlips = max(max(length(WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips), length(WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips))-1, numFlips);
    end
    gapOfInterest = 2;
    
    maxNumCrosses = 0;
    maxTimePoints = 0;
    for flyNum = 1:numFlies
        numCrosses = 0;
        numTimePoints = 0;
        for flipNum = 1:numFlips
            numCrosses = numCrosses + WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).UpProperCrosses(gapOfInterest);
            numCrosses = numCrosses + WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).UpProperCrosses(gapOfInterest);
            numTimePoints = max(length(WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).AbsoluteTime), numTimePoints);
            numTimePoints = max(length(WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).AbsoluteTime), numTimePoints);
        end
        maxNumCrosses = max(numCrosses, maxNumCrosses);
        maxTimePoints = max(numTimePoints, maxTimePoints);
    end
    
    crossBehavMatX = NaN(numFlies, maxNumCrosses, maxTimePoints);
    crossBehavMatY = NaN(numFlies, maxNumCrosses, maxTimePoints);
    for flyNum = 1:numFlies
        crossCounter = 0;
        for flipNum = 1:numFlips
            if WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).UpProperCrosses(gapOfInterest) ~= 0
                for crossInFlip = 1:WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).UpProperCrosses(gapOfInterest)
                    crossCounter = crossCounter + 1;
                    idxForCompID = WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).UpCrossesIndex(gapOfInterest).GapID(crossInFlip);
                    startFrame = WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).UniqCompIDIndex(idxForCompID+1)-1;
                    endFrame = WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).UniqCompIDIndex(idxForCompID+2)-1;
                    crossBehavMatX(flyNum, crossCounter, 1:(endFrame-startFrame+1)) = WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).AlignedFoldedHeadX(startFrame:endFrame) - WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).AlignedFoldedHeadX(startFrame);
                    crossBehavMatY(flyNum, crossCounter, 1:(endFrame-startFrame+1)) = WS.FlipBinnedFlyStruct(flyNum).AlignedData.OddFlips(flipNum).AlignedHeadY(startFrame:endFrame);
                end
            end
            if WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).UpProperCrosses(gapOfInterest) ~= 0
                for crossInFlip = 1:WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).UpProperCrosses(gapOfInterest)
                    crossCounter = crossCounter + 1;
                    idxForCompID = WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).UpCrossesIndex(gapOfInterest).GapID(crossInFlip);
                    startFrame = WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).UniqCompIDIndex(idxForCompID+1)-1;
                    endFrame = WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).UniqCompIDIndex(idxForCompID+2)-1;
                    crossBehavMatX(flyNum, crossCounter, 1:(endFrame-startFrame+1)) = WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).AlignedFoldedHeadX(startFrame:endFrame) - WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).AlignedFoldedHeadX(startFrame);
                    crossBehavMatY(flyNum, crossCounter, 1:(endFrame-startFrame+1)) = WS.FlipBinnedFlyStruct(flyNum).AlignedData.EvenFlips(flipNum).AlignedHeadY(startFrame:endFrame);
                end
            end
        end
        if crossCounter < 5
            crossBehavMatX(flyNum, :, :) = NaN;
            crossBehavMatY(flyNum, :, :) = NaN;
        end
    end

    XYTStruct{WScounter, 1} = crossBehavMatX;
    XYTStruct{WScounter, 2} = crossBehavMatY;

    MeanX = nanmean(max(crossBehavMatX, [], 3), 2);
    StdX = nanstd(max(crossBehavMatX, [], 3), 1, 2);
    avgMeanX = nanmean(nanmean(max(crossBehavMatX, [], 3), 2));
    semMeanX = nanstd(nanmean(max(crossBehavMatX, [], 3), 2))/sqrt(length(nanmean(nanmean(crossBehavMatX, 3), 2)));
    
    % avgMeanXVec(WScounter) = avgMeanX;
    % semMeanXVec(WScounter) = semMeanX;

    MeanXVec(1:numFlies, WScounter) = MeanX;
    StdXVec(1:numFlies, WScounter) = StdX;
    avgMeanXVec(WScounter) = avgMeanX;
    semMeanXVec(WScounter) = semMeanX;
end

% plot(reshape(crossBehavMatX(1,:,:), size(crossBehavMatY,2), size(crossBehavMatY,3)), ...
%      reshape(crossBehavMatY(1,:,:), size(crossBehavMatY,2), size(crossBehavMatY,3)), 'o')
% 
% plot(ones(28,1), nanmean(max(crossBehavMatX, [], 3), 2), 'o');
% plot(ones(28,1), nanmean(nanmean(crossBehavMatX, 3), 2), 'o');

% plot(crossBehavMatX(1,:,:), crossBehavMatY, 'o')

    

% errorbar(avgMeanXVec, semMeanXVec)

labelVec = cell(numWS,1);
for i = 1:numWS
    labelVec{i} = erase(erase(replace(replace(erase(erase(replace(allWS{1, i}.directoryName, '_', ' '), 'lighting '), 'top bottom'), ' IsoD1', ' > +'), 'shts2and3', '> shR'), 'split '), 'new ');
end

figure
hold on
plot(MeanXVec','bo', 'MarkerSize', 3)
plot(nanmean(MeanXVec, 1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
errorbar(1:numWS, nanmean(MeanXVec, 1), nanstd(MeanXVec)/sqrt(sum(~isnan(MeanXVec), 1)-1), 'k', 'CapSize', 0, 'LineWidth', 4)
xlim([0.5, numWS+0.5])
xticklabels(labelVec)
hold off

figure
genoVec = [1, 2];
hold on
plot(MeanXVec(:, genoVec)','bo', 'MarkerSize', 3)
plot(nanmean(MeanXVec(:, genoVec), 1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
errorbar(1:length(genoVec), nanmean(MeanXVec(:, genoVec), 1), nanstd(MeanXVec(:, genoVec))/sqrt(sum(~isnan(MeanXVec(:, genoVec)), 1)-1), 'k', 'CapSize', 0, 'LineWidth', 4)
xlim([0.5, length(genoVec)+.5])
xticks(1:length(genoVec))
xticklabels(labelVec(genoVec))
ylabel('Mean Gap Dip (mm)')
ylim([0, 1])
yticks([0:0.5:1])
title([num2str(gapOfInterest*0.5 + 0.5), ' mm Gap'])
hold off

figure
genoVec = [4, 3, 5];
hold on
plot(MeanXVec(:, genoVec)','bo', 'MarkerSize', 3)
plot(nanmean(MeanXVec(:, genoVec), 1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
errorbar(1:length(genoVec), nanmean(MeanXVec(:, genoVec), 1), nanstd(MeanXVec(:, genoVec))/sqrt(sum(~isnan(MeanXVec(:, genoVec)), 1)-1), 'k', 'CapSize', 0, 'LineWidth', 4)
xlim([0.5, length(genoVec)+.5])
xticks(1:length(genoVec))
xticklabels(labelVec(genoVec))
ylabel('Mean Gap Dip (mm)')
ylim([0, 1])
yticks([0:0.5:1])
title([num2str(gapOfInterest*0.5 + 0.5), ' mm Gap'])
hold off

figure
hold on
for WScounter = [4, 3, 5]
    crossProb = allWS{1, WScounter}.crossStats.AllUpVecProperCrossOverAllGapEventRate(:,gapOfInterest);
    plot(crossProb, MeanXVec(1:length(crossProb), WScounter), 'o')
end
hold off

figure
hold on
genoVec = [1, 2];
pVal = ranksum(mean(sum(~isnan(XYTStruct{genoVec(1),1}), 3),2)/30, mean(sum(~isnan(XYTStruct{genoVec(2),1}), 3),2)/30);
for WScounter = 1:length(genoVec)
    plot(WScounter*ones(size(mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),2))), mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),2)/30, 'bo');
    plot(WScounter, mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),"all")/30, 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
    errorbar(WScounter, mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),"all")/30, std(mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3), 2)/30, 0, 1)/sqrt(size(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),1) - 1), 'k', 'CapSize', 0, 'LineWidth', 4)
    xlim([0.5, length(genoVec)+0.5])
    ylim([0, 0.5])
    yticks([0:0.25:0.5])
    title([num2str(gapOfInterest*0.5 + 0.5), ' mm Gap (p = ', num2str(pVal,5), ')'])
    xticks(1:3)
    xticklabels(labelVec(genoVec))
    ylabel('Gap Dip Duration (s)')
end

figure
hold on
genoVec = [4, 3, 5];
pVal = max(ranksum(mean(sum(~isnan(XYTStruct{genoVec(1),1}), 3),2)/30, mean(sum(~isnan(XYTStruct{genoVec(2),1}), 3),2)/30), ...
           ranksum(mean(sum(~isnan(XYTStruct{genoVec(1),1}), 3),2)/30, mean(sum(~isnan(XYTStruct{genoVec(3),1}), 3),2)/30));
for WScounter = 1:length(genoVec)
    plot(WScounter*ones(size(mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),2))), mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),2)/30, 'bo');
    plot(WScounter, mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),"all")/30, 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
    errorbar(WScounter, mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),"all")/30, std(mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3), 2)/30, 0, 1)/sqrt(size(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),1) - 1), 'k', 'CapSize', 0, 'LineWidth', 4)
    xlim([0.5, length(genoVec)+0.5])
    ylim([0, 0.5])
    yticks([0:0.25:0.5])
    title([num2str(gapOfInterest*0.5 + 0.5), ' mm Gap (p = ', num2str(pVal,5), ')'])
    xticks(1:3)
    xticklabels(labelVec(genoVec))
    ylabel('Gap Dip Duration (s)')
end

figure
hold on
genoVec = [6, 3, 7];
pVal = max(ranksum(mean(sum(~isnan(XYTStruct{genoVec(1),1}), 3),2)/30, mean(sum(~isnan(XYTStruct{genoVec(2),1}), 3),2)/30), ...
           ranksum(mean(sum(~isnan(XYTStruct{genoVec(1),1}), 3),2)/30, mean(sum(~isnan(XYTStruct{genoVec(3),1}), 3),2)/30));
for WScounter = 1:length(genoVec)
    plot(WScounter*ones(size(mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),2))), mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),2)/30, 'bo');
    plot(WScounter, mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),"all")/30, 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
    errorbar(WScounter, mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),"all")/30, std(mean(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3), 2)/30, 0, 1)/sqrt(size(sum(~isnan(XYTStruct{genoVec(WScounter),1}), 3),1) - 1), 'k', 'CapSize', 0, 'LineWidth', 4)
    xlim([0.5, length(genoVec)+0.5])
    ylim([0, 0.5])
    yticks([0:0.25:0.5])
    title([num2str(gapOfInterest*0.5 + 0.5), ' mm Gap (p = ', num2str(pVal,5), ')'])
    xticks(1:3)
    xticklabels(labelVec(genoVec))
    ylabel('Gap Dip Duration (s)')
end

figure
genoVec = [6, 3, 7];
hold on
plot(MeanXVec(:, genoVec)','bo', 'MarkerSize', 3)
plot(nanmean(MeanXVec(:, genoVec), 1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
errorbar(1:3, nanmean(MeanXVec(:, genoVec), 1), nanstd(MeanXVec(:, genoVec))/sqrt(sum(~isnan(MeanXVec(:, genoVec)), 1)-1), 'k', 'CapSize', 0, 'LineWidth', 4)
xlim([0.5, 3.5])
xticks(1:3)
xticklabels(labelVec(genoVec))
ylabel('Mean Gap Dip (mm)')
ylim([0, 1])
yticks([0:0.5:1])
title([num2str(gapOfInterest*0.5 + 0.5), ' mm Gap'])
hold off

figure
hold on
% plot(MeanXVec','bo', 'MarkerSize', 3)
plot(nanmean(MeanXVec./StdXVec, 1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor','k')
% errorbar(1:numWS, nanmean(MeanXVec, 1), nanstd(MeanXVec)/sqrt(sum(~isnan(MeanXVec), 1)-1), 'k', 'CapSize', 0, 'LineWidth', 4)
xlim([0.5, numWS+0.5])
xticklabels(labelVec)
hold off

ranksum(MeanXVec(~isnan(MeanXVec(:,1)),1), MeanXVec(~isnan(MeanXVec(:,2)),2))

% WS.FlipBinnedFlyStruct(3).AlignedData.OddFlips(24).UpCrossesAbsFrameStart(1).GapID(2)