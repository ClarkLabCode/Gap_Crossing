NumFlies = WS.NumFlies;
NumGaps = WS.NumGaps;

% Initialize NumUpCrosses and NumDownCrosses which hold the number of 
% crossing events that happen per fly at each gap width for even and odd 
% flips for both up and down crosses
NumUpCrosses = zeros(NumFlies,2,NumGaps);
NumDownCrosses = zeros(NumFlies,2,NumGaps);

% Fill NumUpCrosses and NumDownCrosses by looping through FBFS and adding 
% the number of crossing events at each gap width per fly for even and odd 
% flips for both up and down crosses
for gapCounter = 1:NumGaps
    for flyCounter = 1:NumFlies
        % Odd flips
        for oddFlipCounter = 1:length(WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips)
            NumUpCrosses(flyCounter,1,gapCounter) = NumUpCrosses(flyCounter,1,gapCounter) + ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrosses(gapCounter);
            NumDownCrosses(flyCounter,1,gapCounter) = NumDownCrosses(flyCounter,1,gapCounter) + ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrosses(gapCounter);
        end
        % Even flips
        for evenFlipCounter = 1:length(WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips)
            NumUpCrosses(flyCounter,2,gapCounter) = NumUpCrosses(flyCounter,2,gapCounter) + ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrosses(gapCounter);
            NumDownCrosses(flyCounter,2,gapCounter) = NumDownCrosses(flyCounter,2,gapCounter) + ...
                WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrosses(gapCounter);
        end
    end
end

% Max number of crossing events at any gap width for any fly of any parity
% This is used as the max size of the fourth dimension of CrossFrameIDMat
MaxNumUpCrosses = max(NumUpCrosses, [], 'all');
MaxNumDownCrosses = max(NumDownCrosses, [], 'all');

% Initialize the matrices that will hold the start and end frame of each
% crossing event for both up and down crosses
% First dimension is Fly#, second is gap# (odd then even), third dimension
% is StartFrame and EndFrame, fourth dimension is crossing event
CrossFrameIDMatUp = zeros(NumFlies,2*NumGaps,2,MaxNumUpCrosses);
CrossFrameIDMatDown = zeros(NumFlies,2*NumGaps,2,MaxNumDownCrosses);

% Now fill CrossFrameIDMat
for gapCounter = 1:NumGaps
    for flyCounter = 1:NumFlies
        % Make a counter that tracks how many crossing events for that
        % given fly at that given gap (first column odd flips, second even)
        CrossCounterUp = ones(NumGaps,2);
        CrossCounterDown = ones(NumGaps,2);
        % Cycle through all odd flips
        for oddFlipCounter = 1:length(WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips)
            % Check to make sure that there are crossing events within the selected flip number
            if WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesIndex(gapCounter).GapID ~= 0
                % Find the first frame of the crossing event
                FirstFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesIndex(gapCounter).GapID+1))-1;
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatUp(flyCounter, gapCounter, 1, (CrossCounterUp(gapCounter,1):(CrossCounterUp(gapCounter,1)+length(FirstFrame)-1))) = ...
                    FirstFrame;
                % Find the last frame of the crossing event
                LastFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesIndex(gapCounter).GapID+2));
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatUp(flyCounter, gapCounter, 2, (CrossCounterUp(gapCounter,1):(CrossCounterUp(gapCounter,1)+length(FirstFrame)-1))) = ...
                    LastFrame;
                % Advance the cross counter appropriately
                CrossCounterUp(gapCounter,1) = CrossCounterUp(gapCounter,1)+length(FirstFrame);
            end
            % Check to make sure that there are crossing events within the selected flip number
            if WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesIndex(gapCounter).GapID ~= 0
                % Find the first frame of the crossing event
                FirstFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesIndex(gapCounter).GapID+1))-1;
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatDown(flyCounter, gapCounter, 1, (CrossCounterDown(gapCounter,1):(CrossCounterDown(gapCounter,1)+length(FirstFrame)-1))) = ...
                    FirstFrame;
                % Find the last frame of the crossing event
                LastFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesIndex(gapCounter).GapID+2));
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatDown(flyCounter, gapCounter, 2, (CrossCounterDown(gapCounter,1):(CrossCounterDown(gapCounter,1)+length(FirstFrame)-1))) = ...
                    LastFrame;
                % Advance the cross counter appropriately
                CrossCounterDown(gapCounter,1) = CrossCounterDown(gapCounter,1)+length(FirstFrame);
            end
        end
        
        % Cycle through all even flips
        for evenFlipCounter = 1:length(WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips)
            % Check to make sure that there are crossing events within the selected flip number
            if WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesIndex(gapCounter).GapID ~= 0
                % Find the first frame of the crossing event
                FirstFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesIndex(gapCounter).GapID+1))-1;
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatUp(flyCounter, (gapCounter + NumGaps), 1, (CrossCounterUp(gapCounter,2):(CrossCounterUp(gapCounter,2)+length(FirstFrame)-1))) = ...
                    FirstFrame;
                % Find the last frame of the crossing event
                LastFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesIndex(gapCounter).GapID+2));
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatUp(flyCounter, (gapCounter + NumGaps), 2, (CrossCounterUp(gapCounter,2):(CrossCounterUp(gapCounter,2)+length(FirstFrame)-1))) = ...
                    LastFrame;
                % Advance the cross counter appropriately
                CrossCounterUp(gapCounter,2) = CrossCounterUp(gapCounter,2)+length(FirstFrame);
            end
            % Check to make sure that there are crossing events within the selected flip number
            if WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesIndex(gapCounter).GapID ~= 0
                % Find the first frame of the crossing event
                FirstFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesIndex(gapCounter).GapID+1))-1;
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatDown(flyCounter, (gapCounter + NumGaps), 1, (CrossCounterDown(gapCounter,2):(CrossCounterDown(gapCounter,2)+length(FirstFrame)-1))) = ...
                    FirstFrame;
                % Find the last frame of the crossing event
                LastFrame = ...
                    WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).AbsoluteTime(...
                        WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).UniqCompIDIndex(...
                            WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesIndex(gapCounter).GapID+2));
                % Now fill in the corresponding value in CrossFrameIDMat
                % If there were multiple events at a given width for a
                % given fly within a single flip, the : vector notation
                % fits it in appropriately
                CrossFrameIDMatDown(flyCounter, (gapCounter + NumGaps), 2, (CrossCounterDown(gapCounter,2):(CrossCounterDown(gapCounter,2)+length(FirstFrame)-1))) = ...
                    LastFrame;
                % Advance the cross counter appropriately
                CrossCounterDown(gapCounter,2) = CrossCounterDown(gapCounter,2)+length(FirstFrame);
            end
        end
    end
end