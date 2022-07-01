% Takes the probabilities output by NN for "crossing" events and uses them
% to classify each event as either a Glass Crossing or Proper Crossing

% This function has a built in threshold that is used for classification.
% This threshold should consistently remain the same between analyses
% unless showing that the threshold has little-to-no effect on the results.
% The default for this threshold value should be 0.5

% This function also renormalizes the probabilities such that:
% P(On)     = P(On)     /   (P(On) + P(Off))
% P(Off)    = P(Off)    /   (P(On) + P(Off))
% This effectively sets P(Ambig) = 0 and redistributes that probability
% This is done because this was found to produce the best NN classification
% results during the NN training procedures done by Joseph and Baohua

function WS = ClassifyOnOff(WS)

% Port in the relevant fields from WS
NumGaps = WS.NumGaps;
threshProb = WS.threshProb;

% Produce the renormalized probabilities such that:
% P(On)     = P(On)     /   (P(On) + P(Off))
% P(Off)    = P(Off)    /   (P(On) + P(Off))
for flyCounter = 1:length(WS.FlipBinnedFlyStruct.ExpNum)
    % Do this for odd flips
    for flipCounter = 1:length(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips)
        for gapCounter = 1:NumGaps
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).RenormUpGlassCrossProb(gapCounter).GapID = ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesOnProb(gapCounter).GapID ./ ...
            (WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesOnProb(gapCounter).GapID + ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrossesOffProb(gapCounter).GapID);
        
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).RenormDownGlassCrossProb(gapCounter).GapID = ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesOnProb(gapCounter).GapID ./ ...
            (WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesOnProb(gapCounter).GapID + ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrossesOffProb(gapCounter).GapID);
        end
    end
    % Do this for even flips
    for flipCounter = 1:length(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips)
        for gapCounter = 1:NumGaps
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).RenormUpGlassCrossProb(gapCounter).GapID = ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesOnProb(gapCounter).GapID ./ ...
            (WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesOnProb(gapCounter).GapID + ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrossesOffProb(gapCounter).GapID);
        
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).RenormDownGlassCrossProb(gapCounter).GapID = ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesOnProb(gapCounter).GapID ./ ...
            (WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesOnProb(gapCounter).GapID + ...
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrossesOffProb(gapCounter).GapID);
        end
    end
end

% Now go through with the classifications of each "crossing" event
% Do this by checking if the threshold probability is met or not
% The sum( ) is there because there can be multiple "crossing" events at a
% given gap width within one flip and we care for the total number of them
for flyCounter = 1:length(WS.FlipBinnedFlyStruct.ExpNum)
    % Odd flips
    for flipCounter = 1:length(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips)
        for gapCounter = 1:NumGaps
            % Do this for up and down crossings
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpProperCrosses(gapCounter) = ...
                sum(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).RenormUpGlassCrossProb(gapCounter).GapID<threshProb);
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpGlassCrosses(gapCounter) = ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpCrosses(gapCounter) - ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).UpProperCrosses(gapCounter);
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownProperCrosses(gapCounter) = ...
                sum(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).RenormDownGlassCrossProb(gapCounter).GapID<threshProb);
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownGlassCrosses(gapCounter) = ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownCrosses(gapCounter) - ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(flipCounter).DownProperCrosses(gapCounter);
        end
    end
    % Even flips
    for flipCounter = 1:length(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips)
        for gapCounter = 1:NumGaps
            % Do this for up and down crossings
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpProperCrosses(gapCounter) = ...
                sum(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).RenormUpGlassCrossProb(gapCounter).GapID<threshProb);
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpGlassCrosses(gapCounter) = ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpCrosses(gapCounter) - ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).UpProperCrosses(gapCounter);
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownProperCrosses(gapCounter) = ...
                sum(WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).RenormDownGlassCrossProb(gapCounter).GapID<threshProb);
            WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownGlassCrosses(gapCounter) = ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownCrosses(gapCounter) - ...
                WS.FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips(flipCounter).DownProperCrosses(gapCounter);
        end
    end
end

WS.NNClassificationComplete = 'Yes';

end
