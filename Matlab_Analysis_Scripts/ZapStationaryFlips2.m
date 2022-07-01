% Deletes data within flips where flies spend any time fully stationary

function WS = ZapStationaryFlips2(WS)

% function FlipBinnedFlyStruct = ZapStationaryFlips(FlipBinnedFlyStruct)

% Port in the relevant fields from WS
FlipBinnedFlyStruct = WS.FlipBinnedFlyStruct;
PreZapFBFS = FlipBinnedFlyStruct;

% Initialize variable that holds the max length of flips in frames
maxFramesPerFlip = 0;

% Find the max number of frames in a flip by finding the max
% number of frames for each fly and taking the max of that
for expCounter = 1:length(FlipBinnedFlyStruct)
    for flyCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum)
        
        % Find the max number of frames in an odd flip and store it
        for oddFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips)
            maxFramesPerFlip = max(length(...
                [FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime]),...
                maxFramesPerFlip);
        end
        
        % Find the max number of frames in an even flip and store it
        for evenFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips)
            maxFramesPerFlip = max(length(...
                [FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime]),...
                maxFramesPerFlip);
        end
    end
end
        
% Cycle through every fly in every experiment to check where the fly is too
% inactive in a flip and set it to 0 or where the fly is lost track of
% during a flip for moer than 5 consecutive frames
for expCounter = 1:length(FlipBinnedFlyStruct)
    for flyCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum)
        for oddFlipCounter = 1:length(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips)
            % If there are flips in which more than 5 consecutive frames are lost or in 
            % which the fly is inactive for at least half the flip, replace the data with 0
            if (max([diff(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime),0]) > 5) ...
                    || length([FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime]) ...
                        < (maxFramesPerFlip/2)
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
            % If there are flips in which more than 5 consecutive frames are lost or in 
            % which the fly is inactive for at least half the flip, replace the data with 0
            if (max([diff(FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime),0]) > 5) ...
                    || length([FlipBinnedFlyStruct(expCounter).ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime]) ...
                        < (maxFramesPerFlip/2)
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

% Update the fields in WS
WS.FlipBinnedFlyStruct = FlipBinnedFlyStruct;
WS.PreZapFBFS = PreZapFBFS;

end