NumCorridors = 7;
NumGaps = 4;

% Initialize NumCrosses which holds the number of crossing events that
% happen per fly at each gap width for even and odd flips
NumCrosses = zeros(NumCorridors,2,NumGaps);

% Fill NumCrosses by looping through LabeledFlyStruct and adding the number
% of crossing events at ech gap width per fly for even and odd flips
for gapCounter = 1:NumGaps
    for flyCounter = 1:NumCorridors
        % Odd flips
        for oddFlipCounter = 1:length(LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips)
            NumCrosses(flyCounter,1,gapCounter) = NumCrosses(flyCounter,1,gapCounter) + ...
                LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).Crosses(gapCounter);
        end
        % Even flips
        for evenFlipCounter = 1:length(LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips)
            NumCrosses(flyCounter,2,gapCounter) = NumCrosses(flyCounter,2,gapCounter) + ...
                LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).Crosses(gapCounter);
        end
    end
end

% Max number of crossing events at any gap width for any fly of any parity
% This is used as the max size of the fourth dimension of CrossFrameIDMat
MaxNumCrosses = max(NumCrosses, [], 'all');

% Initialize the matrix that will hold the start and end frame of each
% crossing event
% First dimension is Fly#, second is gap# (odd then even), third dimension
% is StartFrame and EndFrame, fourth dimension is crossing event
CrossFrameIDMat = zeros(NumCorridors,2*NumGaps,2,MaxNumCrosses);

% Now fill CrossFrameIDMat
for gapCounter = 1:NumGaps
    for flyCounter = 1:NumCorridors
        % Make a counter that tracks how many crossing events for that
        % given fly at that given gap (first column odd flips, second even)
        CrossCounter = ones(NumGaps,2);
        % Cycle through all odd flips
        for oddFlipCounter = 1:length(LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips)
            % Check to make sure that there are crossing events within the selected flip number
            if LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).CrossesIndex(gapCounter).GapID ~= 0
                % Find the first frame of the crossing event
                FirstFrame = ...
                    LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime(...
                        LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UniqCompIDIndex(...
                            LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).CrossesIndex(gapCounter).GapID+1))-1;
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMat(flyCounter, gapCounter, 1, (CrossCounter(gapCounter,1):(CrossCounter(gapCounter,1)+length(FirstFrame)-1))) = ...
                    FirstFrame;
                % Find the last frame of the crossing event
                LastFrame = ...
                    LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime(...
                        LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UniqCompIDIndex(...
                            LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).CrossesIndex(gapCounter).GapID+2));
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMat(flyCounter, gapCounter, 2, (CrossCounter(gapCounter,1):(CrossCounter(gapCounter,1)+length(FirstFrame)-1))) = ...
                    LastFrame;
                % Advance the cross counter appropriately
                CrossCounter(gapCounter,1) = CrossCounter(gapCounter,1)+length(FirstFrame);
            end
        end
        
        % Cycle through all even flips
        for evenFlipCounter = 1:length(LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips)
            % Check to make sure that there are crossing events within the selected flip number
            if LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).CrossesIndex(gapCounter).GapID ~= 0
                % Find the first frame of the crossing event
                FirstFrame = ...
                    LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime(...
                        LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UniqCompIDIndex(...
                            LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).CrossesIndex(gapCounter).GapID+1))-1;
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMat(flyCounter, (gapCounter + NumGaps), 1, (CrossCounter(gapCounter,2):(CrossCounter(gapCounter,2)+length(FirstFrame)-1))) = ...
                    FirstFrame;
                % Find the last frame of the crossing event
                LastFrame = ...
                    LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime(...
                        LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UniqCompIDIndex(...
                            LabeledFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).CrossesIndex(gapCounter).GapID+2));
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMat(flyCounter, (gapCounter + NumGaps), 2, (CrossCounter(gapCounter,2):(CrossCounter(gapCounter,2)+length(FirstFrame)-1))) = ...
                    LastFrame;
                % Advance the cross counter appropriately
                CrossCounter(gapCounter,2) = CrossCounter(gapCounter,2)+length(FirstFrame);
            end
        end
    end
end