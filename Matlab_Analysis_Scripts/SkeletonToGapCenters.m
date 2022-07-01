% Uses the user-made skeleton to find the centers of the gaps

function WS = SkeletonToGapCenters(WS)

% Port in the relevant fields from WS
NumGaps = WS.NumGaps;
NumCorridors = WS.NumCorridors;
Skeleton_x = WS.Skeleton_x;
Skeleton_y = WS.Skeleton_y;

% Initialize the matices that hold the positions of the gap centers
GapCentsXOdd  = zeros(NumGaps,NumCorridors);
GapCentsYOdd  = zeros(NumGaps,NumCorridors);
GapCentsXEven = zeros(NumGaps,NumCorridors);
GapCentsYEven = zeros(NumGaps,NumCorridors);

% Loop through all corridors and gaps to find the centers in each orientation
% Find the centers by averaging two diagonally opposite points in the gaps
for corrCounter = 1:NumCorridors
    for gapCounter = 1:NumGaps
        GapCentsXOdd(gapCounter,corrCounter) = ...
            mean([Skeleton_x(4*(gapCounter-1)+2,1,corrCounter),...
                  Skeleton_x(4*(2*NumGaps-gapCounter)+4,1,corrCounter)]);
        GapCentsYOdd(gapCounter,corrCounter) = ...
            mean([Skeleton_y(4*(gapCounter-1)+2,1,corrCounter),...
                  Skeleton_y(4*(2*NumGaps-gapCounter)+4,1,corrCounter)]);
        GapCentsXEven(gapCounter,corrCounter) = ...
            mean([Skeleton_x(4*(gapCounter-1)+2,2,corrCounter),...
                  Skeleton_x(4*(2*NumGaps-gapCounter)+4,2,corrCounter)]);
        GapCentsYEven(gapCounter,corrCounter) = ...
            mean([Skeleton_y(4*(gapCounter-1)+2,2,corrCounter),...
                  Skeleton_y(4*(2*NumGaps-gapCounter)+4,2,corrCounter)]);
    end
end

% Round the positions of the gap centers so that it's an integer pixel
GapCentsXOdd    = round(GapCentsXOdd);
GapCentsYOdd    = round(GapCentsYOdd);
GapCentsXEven   = round(GapCentsXEven);
GapCentsYEven   = round(GapCentsYEven);

% Update the fields in WS
WS.GapCentsXOdd = GapCentsXOdd;
WS.GapCentsYOdd = GapCentsYOdd;
WS.GapCentsXEven = GapCentsXEven;
WS.GapCentsYEven = GapCentsYEven;

end