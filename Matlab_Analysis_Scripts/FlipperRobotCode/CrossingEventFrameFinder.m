% Adds fields to FBFS that contain the start and end frame of each
% "crossing" event per fly at all flips for each gap width

function WS = CrossingEventFrameFinder(WS)

% Port in the relevant fields from WS
FBFS = WS.FlipBinnedFlyStruct;
NumGaps = WS.NumGaps;

% Below is how to get the index of the crossing event within up crosses
% FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipNum).UpCrossesIndex(gapNum).GapID

% Feed this into below to get the frameNum within flip of the crossing event
% FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipNum).UniqCompIDIndex(__)

% Now feed that into below to get the absolute frameNum of the event
% FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipNum).AbsoluteTime(__)

% The above gets the first frame at which the fly enters the crossing
% sequence, but this isn't the first frame of the crossing event

% First frame of the crossing event happens when you add 1 to the index
% (first part you retrieve above) and then subtract 1 from the absolute time
% (third part you retrieve above)

% This is because the event is when a transition such as 3->2->1 occurs,
% but the crossing is the first frame of being in 2 until last frame of
% being in 2

% Note that we need to do this for up/down and even/odd and also must loop
% through all flies, all flips, and all gaps

% Do all odd flips first in this nest of loops
for flyCounter = 1:length(FBFS.ExpNum)
    for flipCounter = 1:length(FBFS.ExpNum(flyCounter).BehavData.OddFlips)
        for gapCounter = 1:NumGaps
            % Odd and Up start frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesAbsFrameStart(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesIndex(gapCounter).GapID + 1 ...
                                                                                             ) ...
                                                                                      ) - 1;
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesAbsFrameStart(gapCounter).GapID = 0;
            end
            % Odd and Down start frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesAbsFrameStart(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesIndex(gapCounter).GapID + 1 ...
                                                                                             ) ...
                                                                                      ) - 1; 
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesAbsFrameStart(gapCounter).GapID = 0;
            end
            % Odd and Up end frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesAbsFrameEnd(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesIndex(gapCounter).GapID + 2 ...
                                                                                             ) ...
                                                                                      ) + 1;
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesAbsFrameEnd(gapCounter).GapID = 0;
            end
            % Odd and Down end frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesAbsFrameEnd(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesIndex(gapCounter).GapID + 2 ...
                                                                                             ) ...
                                                                                      ) + 1;  
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesAbsFrameEnd(gapCounter).GapID = 0;
            end
            
        end
    end
end

% Now do all even flips in this nest of loops
for flyCounter = 1:length(FBFS.ExpNum)
    for flipCounter = 1:length(FBFS.ExpNum(flyCounter).BehavData.EvenFlips)
        for gapCounter = 1:NumGaps
            % Even and Up start frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesAbsFrameStart(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesIndex(gapCounter).GapID + 1 ...
                                                                                             ) ...
                                                                                      ) - 1;    
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesAbsFrameStart(gapCounter).GapID = 0;
            end
            % Even and Down start frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesAbsFrameStart(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesIndex(gapCounter).GapID + 1 ...
                                                                                             ) ...
                                                                                      ) - 1;  
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesAbsFrameStart(gapCounter).GapID = 0;
            end
            % Even and Up end frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesAbsFrameEnd(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesIndex(gapCounter).GapID + 2 ...
                                                                                             ) ...
                                                                                      ) + 1;  
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesAbsFrameEnd(gapCounter).GapID = 0;
            end
            % Even and Down end frame
            % If there's a crossing
            if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesIndex(gapCounter).GapID
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesAbsFrameEnd(gapCounter).GapID = ...
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).AbsoluteTime( ...
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UniqCompIDIndex( ...
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesIndex(gapCounter).GapID + 2 ...
                                                                                             ) ...
                                                                                      ) + 1; 
            % If there isn't a crossing
            else
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesAbsFrameEnd(gapCounter).GapID = 0;
            end
                                                                  
        end
    end
end

% Update the fields in WS
WS.FlipBinnedFlyStruct = FBFS;

end
