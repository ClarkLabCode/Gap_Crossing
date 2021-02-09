% Asks the user to outline where each corridor is within the video, one
% time for each orientation of the cassette
% Then identifies which corridor each row of finalStats is in and appends
% it to the structure

% Inputs
% inputFileName     = Name of video file to be analyzed, must include extension
% finalStats    = Structure that holds all info from the video, each row
%                 corresponds to a particular fly at a particular time
% NumCorridors  = Number of corridors in video

% Outputs
% finalStats    = Structure that holds all info from the video, each row
%                 corresponds to a particular fly at a particular time
% CorrMask1     = Cell array that holds the mask for each corridor for all
%                 odd numbered flips
% CorrMask2     = Cell array that holds the mask for each corridor for all
%                 even numbered flips

function [finalStats, CorrMask1, CorrMask2] = ...
    CorridorIdentifier(inputFileName, finalStats, NumCorridors)

% Declare VideoReader object
reader1 = VideoReader(inputFileName);

% CorrMask holds within it the mask for each corridor, and 1/2 correspond
% to the orientation of the cassettes, with 1 being the odd numbered flips
% and 2 being the even numbered flips
CorrMask1 = cell(1, NumCorridors);
CorrMask2 = cell(1, NumCorridors);

% Let the user draw the mask for each corridor and save in the appropriate
% element within CorrMask1/2. Gives text output to ensure user keeps track
% of which corridor they must label.
for CorrCounter = 1:NumCorridors
    fprintf(strcat('Select Corridor #',num2str(CorrCounter),'in orientation 1\n'));
    CorrMask1{CorrCounter} = roipoly(read(reader1,indPos(1)+5));
    fprintf(strcat('Select Corridor #',num2str(CorrCounter),'in orientation 2\n'));
    CorrMask2{CorrCounter} = roipoly(read(reader1,indPos(3)+5)); 
end

% Go through finalStats and determine which corridor each fly is in
for Row = 1:length(finalStats)
    CentroidX = round(finalStats(Row).Centroid(1));
    CentroidY = round(finalStats(Row).Centroid(2));
    finalStats(Row).CorridorID = 0;     % Catches empty cells
    
    % Odd numbered flips
    if ((-1)^(finalStats(Row).FlipNumber) ~= 1)
        for i = 1:NumCorridors
            if CorrMask1{i}(CentroidY, CentroidX) == 1
                finalStats(Row).CorridorID = i;
                break
            end
        end
    % Even numbered flips
    else 
        for i = 1:NumCorridors
            if CorrMask2{i}(CentroidY, CentroidX) == 1
                finalStats(Row).CorridorID = i;
                break
            end
        end 
    end
end

end