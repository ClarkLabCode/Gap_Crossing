% Odd flips
cd ..\..\Data\IsoD1_test\Background_Frames
figure(2)
imshow(imread('IsoD1_test_bg_flip1.png'));
cd ..\..\..\Matlab_Analysis_Scripts\FlipperRobotCode

hold on

flyCounter = 1;

for flipCounter = 1:length([FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips])
    for flyCounter = 1:length([FlipBinnedFlyStruct.ExpNum])
        if ~isempty(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).CentroidX)
            plot(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).CentroidX(1),...
                FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).CentroidY(1), 'bo','LineWidth',1);
            plot(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).CentroidX,...
                FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).CentroidY, 'r');
        end
    end
end    

hold off

% Even flips
cd ..\..\Data\IsoD1_test\Background_Frames
figure(1)
imshow(imread('IsoD1_test_bg_flip2.png'));
cd ..\..\..\Matlab_Analysis_Scripts\FlipperRobotCode

hold on

for flipCounter = 1:length([FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips])
    for flyCounter = 1:length([FlipBinnedFlyStruct.ExpNum])
        if ~isempty(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).CentroidX)
            plot(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).CentroidX(1),...
                FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).CentroidY(1), 'bo','LineWidth',1);
            plot(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).CentroidX,...
                FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).CentroidY, 'r');
        end
    end
end

hold off