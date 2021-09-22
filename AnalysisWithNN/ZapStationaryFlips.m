% Deletes data within flips where flies spend any time fully stationary

function FlipBinnedFlyStruct = ZapStationaryFlips(FlipBinnedFlyStruct)

for expCounter = 1:length(FlipBinnedFlyStruct)
    for flyCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum)
        maxFramesPerOddFlip = 0;
        maxFramesPerEvenFlip = 0;
        for oddFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips)
            % Checks for flip numbers in which there are fewer frames than the maximum expected amount
            maxFramesPerOddFlip = max(length(...
                [FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime]),...
                maxFramesPerOddFlip);
        end
        for evenFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips)
            % Checks for flip numbers in which there are fewer frames than the maximum expected amount
            maxFramesPerEvenFlip = max(length(...
                [FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime]),...
                maxFramesPerEvenFlip);
        end
        for oddFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips)
            % If there are fewer frames in a flip than expected, replace the data with 0
            if length([FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime]) ...
                    < (maxFramesPerOddFlip-1)
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).Area = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).CentroidX = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).CentroidY = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).FrameInFlip = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).BoundingBoxTLX = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).BoundingBoxTLY = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).BoundingBoxWidth = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).BoundingBoxHeight = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).MajorAxisLength = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).MinorAxisLength = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).Circularity = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).Orientation = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).VelX = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).VelY = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).Speed = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AngularSpeed = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).CompID = 0;
            end
        end
        for evenFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips)
            % If there are fewer frames in a flip than expected, replace the data with 0
            if length([FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime]) ...
                    < (maxFramesPerEvenFlip-1)
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).Area = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).CentroidX = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).CentroidY = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).FrameInFlip = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).BoundingBoxTLX = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).BoundingBoxTLY = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).BoundingBoxWidth = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).BoundingBoxHeight = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).MajorAxisLength = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).MinorAxisLength = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).Circularity = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).Orientation = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).VelX = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).VelY = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).Speed = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AngularSpeed = 0;
                FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).CompID = 0;
            end
        end
    end
end