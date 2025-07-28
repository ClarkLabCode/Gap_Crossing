fileName = 'C:\Users\clark\Joe\Gap_Crossing\Data\IsoD1_lighting_top_bottom\Experiment_1\0_Raw_Videos\2022-06-26 11-51-52_IsoD1_lighting_top_bottom_Exp1.mov';
vR = VideoReader(fileName);
vW = VideoWriter('Supp_Movie.avi');
vW.FrameRate = 30;


timeVec = [];
xVec = [];
yVec = [];

for flipCounter = 4:6
    for flyCounter = 1:7
        timeVec = [timeVec, WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(flipCounter).AbsoluteTime];
        timeVec = [timeVec, WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(flipCounter).AbsoluteTime];
        xVec = [xVec, WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(flipCounter).CentroidX];
        xVec = [xVec, WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(flipCounter).CentroidX];
        yVec = [yVec, WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(flipCounter).CentroidY];
        yVec = [yVec, WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(flipCounter).CentroidY];
    end
end

[sortedTimeVec, sortingIdx] = sort(timeVec);
sortingIdx = sortingIdx(sortedTimeVec~=0);
sortedTimeVec = timeVec(sortingIdx);
sortedXVec = xVec(sortingIdx);
sortedYVec = yVec(sortingIdx);

sum(diff(sortedTimeVec)>2);
idxOfFlips = [1, find(diff(sortedTimeVec)>2), length(sortedTimeVec)];
allTimes1 = unique(sortedTimeVec(idxOfFlips(1):idxOfFlips(2)));
allTimes2 = unique(sortedTimeVec(idxOfFlips(2)+1:idxOfFlips(3)));
allTimes3 = unique(sortedTimeVec(idxOfFlips(3)+1:idxOfFlips(4)));
allTimes4 = unique(sortedTimeVec(idxOfFlips(4)+1:idxOfFlips(5)));
allTimes5 = unique(sortedTimeVec(idxOfFlips(5)+1:idxOfFlips(6)));
allTimes6 = unique(sortedTimeVec(idxOfFlips(6)+1:idxOfFlips(7)));
allTimes = unique(sortedTimeVec);

plotTime1 = NaN(7, length(allTimes1));
plotX1 = NaN(7, length(allTimes1));
plotY1 = NaN(7, length(allTimes1));

plotTime2 = NaN(7, length(allTimes2));
plotX2 = NaN(7, length(allTimes2));
plotY2 = NaN(7, length(allTimes2));

plotTime3 = NaN(7, length(allTimes3));
plotX3 = NaN(7, length(allTimes3));
plotY3 = NaN(7, length(allTimes3));

plotTime4 = NaN(7, length(allTimes4));
plotX4 = NaN(7, length(allTimes4));
plotY4 = NaN(7, length(allTimes4));

for flyCounter = 1:7
    for timeCounter = 1:length(allTimes1)
        if sum(WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).AbsoluteTime == allTimes1(timeCounter)) ~= 0
            plotTime1(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).AbsoluteTime(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).AbsoluteTime == allTimes1(timeCounter));
            plotX1(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).CentroidX(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).AbsoluteTime == allTimes1(timeCounter));
            plotY1(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).CentroidY(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(4).AbsoluteTime == allTimes1(timeCounter));
        end
    end
end

for flyCounter = 1:7
    for timeCounter = 1:length(allTimes2)
        if sum(WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).AbsoluteTime == allTimes2(timeCounter)) ~= 0
            plotTime2(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).AbsoluteTime(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).AbsoluteTime == allTimes2(timeCounter));
            plotX2(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).CentroidX(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).AbsoluteTime == allTimes2(timeCounter));
            plotY2(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).CentroidY(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(4).AbsoluteTime == allTimes2(timeCounter));
        end
    end
end

for flyCounter = 1:7
    for timeCounter = 1:length(allTimes3)
        if sum(WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).AbsoluteTime == allTimes3(timeCounter)) ~= 0
            plotTime3(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).AbsoluteTime(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).AbsoluteTime == allTimes3(timeCounter));
            plotX3(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).CentroidX(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).AbsoluteTime == allTimes3(timeCounter));
            plotY3(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).CentroidY(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(5).AbsoluteTime == allTimes3(timeCounter));
        end
    end
end

for flyCounter = 1:7
    for timeCounter = 1:length(allTimes4)
        if sum(WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).AbsoluteTime == allTimes4(timeCounter)) ~= 0
            plotTime4(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).AbsoluteTime(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).AbsoluteTime == allTimes4(timeCounter));
            plotX4(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).CentroidX(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).AbsoluteTime == allTimes4(timeCounter));
            plotY4(flyCounter,timeCounter) = ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).CentroidY(...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(5).AbsoluteTime == allTimes4(timeCounter));
        end
    end
end

open(vW);
for timeCounter = 1:length(allTimes1)
    frameImg = read(vR, allTimes1(timeCounter));
    imshow(medfilt2(frameImg(:,:,1)));
    hold on
    for flyCounter = 1:7
        plot(plotX1(flyCounter,timeCounter), plotY1(flyCounter,timeCounter), 'r*')
        plot(plotX1(flyCounter,1:timeCounter), plotY1(flyCounter,1:timeCounter), 'r')
    end
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end

for timeCounter = (allTimes1(end)+1):(allTimes2(1)-1)
    frameImg = read(vR, timeCounter);
    imshow(medfilt2(frameImg(:,:,1)));
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end

for timeCounter = 1:length(allTimes2)
    frameImg = read(vR, allTimes2(timeCounter));
    imshow(medfilt2(frameImg(:,:,1)));
    hold on
    for flyCounter = 1:7
        plot(plotX2(flyCounter,timeCounter), plotY2(flyCounter,timeCounter), 'r*')
        plot(plotX2(flyCounter,1:timeCounter), plotY2(flyCounter,1:timeCounter), 'r')
    end
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end

for timeCounter = (allTimes2(end)+1):(allTimes3(1)-1)
    frameImg = read(vR, timeCounter);
    imshow(medfilt2(frameImg(:,:,1)));
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end

for timeCounter = 1:length(allTimes3)
    frameImg = read(vR, allTimes3(timeCounter));
    imshow(medfilt2(frameImg(:,:,1)));
    hold on
    for flyCounter = 1:7
        plot(plotX3(flyCounter,timeCounter), plotY3(flyCounter,timeCounter), 'r*')
        plot(plotX3(flyCounter,1:timeCounter), plotY3(flyCounter,1:timeCounter), 'r')
    end
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end

for timeCounter = (allTimes3(end)+1):(allTimes4(1)-1)
    frameImg = read(vR, timeCounter);
    imshow(medfilt2(frameImg(:,:,1)));
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end

for timeCounter = 1:length(allTimes4)
    frameImg = read(vR, allTimes4(timeCounter));
    imshow(medfilt2(frameImg(:,:,1)));
    hold on
    for flyCounter = 1:7
        plot(plotX4(flyCounter,timeCounter), plotY4(flyCounter,timeCounter), 'r*')
        plot(plotX4(flyCounter,1:timeCounter), plotY4(flyCounter,1:timeCounter), 'r')
    end
    drawnow()
    frameToSave = frame2im(getframe(gcf));
    croppedFrame = frameToSave(100:713,100:1354,:);
    writeVideo(vW,croppedFrame);
end
close(vW);