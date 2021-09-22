% Figures out which compartment each row of finalStats is in and appends it
% to the structure

% Inputs:
% finalStats    = Structure that holds all info from the video, each row
%                 corresponds to a particular fly at a particular time
% CompMask1     = Cell array that holds the mask for each compartment for all
%                 odd numbered flips
% CompMask2     = Cell array that holds the mask for each compartment for all
%                 even numbered flips
% indPos        = Index of frame at which flips happen
% NumComps      = Number of compartments per corridor

% Output:
% finalStats    = Structure that holds all info from the video, each row
%                 corresponds to a particular fly at a particular time

function finalStats = CompIdentifier(finalStats, CompMask1, CompMask2, indPos, NumComps)

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

end