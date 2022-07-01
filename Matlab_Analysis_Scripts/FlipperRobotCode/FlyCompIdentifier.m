% Figures out which compartment each row of finalStats is in and appends it
% to the structure

function WS = FlyCompIdentifier(WS)

% function finalStats = CompIdentifier(finalStats, CompMask1, CompMask2, indPos, NumComps)

% Port in the relevant fields from WS
finalStats          = WS.finalStats;
CompMask1           = WS.CompMask1;
CompMask2           = WS.CompMask2;
indPos              = WS.indPos;
NumComps            = WS.NumComps;

% Go through each row of finalStats and figure out which compartment the
% fly is in and add that info to the structure
for Row = 1:length(finalStats)
    CentroidX = round(finalStats(Row).Centroid(1));
    CentroidY = round(finalStats(Row).Centroid(2));
    finalStats(Row).CompID = 0;     % Catches empty cells
    % Even flip numbers
    if ((-1)^(finalStats(Row).FlipNumber) ~= 1)
        for i = 1:NumComps
            if CompMask1{i}(CentroidY, CentroidX) == 1
                finalStats(Row).CompID = i;
                break
            end
        end
    % Odd flip numbers
    else
        for i = 1:NumComps
            if CompMask2{i}(CentroidY, CentroidX) == 1
                finalStats(Row).CompID = i;
                break
            end
        end 
    end
end

% Update the fields in WS
WS.finalStats = finalStats;

end
