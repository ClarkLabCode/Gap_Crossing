% Function that uses the previously found corridor masks to identify which
% corridor each fly is in inside of the finalStats struct

function WS = FlyCorrIdentifier(WS)

% Port in the relevant fields from WS
finalStats          = WS.finalStats;
NumCorridors        = WS.NumCorridors;
CorrMask1           = WS.CorrMask1;
CorrMask2           = WS.CorrMask2;

% Go through finalStats and determine which corridor each fly is in
for Row = 1:length(finalStats)
    CentroidX = round(finalStats(Row).Centroid(1));
    CentroidY = round(finalStats(Row).Centroid(2));
    finalStats(Row).CorridorID = 0;     % Catches empty cells

    % Odd numbered flips
    if ((-1)^(finalStats(Row).FlipNumber) ~= 1)
        for corrCounter = 1:NumCorridors
            if CorrMask1{corrCounter}(CentroidY, CentroidX) == 1
                finalStats(Row).CorridorID = corrCounter;
                break
            end
        end
    % Even numbered flips
    else 
        for corrCounter = 1:NumCorridors
            if CorrMask2{corrCounter}(CentroidY, CentroidX) == 1
                finalStats(Row).CorridorID = corrCounter;
                break
            end
        end 
    end
end

% Update the fields in WS
WS.finalStats = finalStats;

end