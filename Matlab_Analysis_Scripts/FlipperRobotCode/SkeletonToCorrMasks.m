% Uses the user-made skeleton to create the mask for the corridors

function WS = SkeletonToCorrMasks(WS)

% Port in the relevant fields from WS
NumGaps = WS.NumGaps;
NumCorridors = WS.NumCorridors;
Skeleton_x = WS.Skeleton_x;
Skeleton_y = WS.Skeleton_y;

% Initialize the mask cells
CorrMask1 = cell(1, NumCorridors);
CorrMask2 = cell(1, NumCorridors);

% Initialize a temporary mask to be used
tempMask1 = false(1080,1920);
tempMask2 = false(1080,1920);

% Loop through the skeleton data for each corridor and fill in for each orientation
% For poly2mask to work, we need to close the polygon by repeating the first point
% We then dilate the mask slightly to make sure we capture the full corridor
for corrCounter = 1:NumCorridors
    tempMask1 = poly2mask([Skeleton_x(:,1,corrCounter); Skeleton_x(1,1,corrCounter)],...
        [Skeleton_y(:,1,corrCounter);Skeleton_y(1,1,corrCounter)],1080,1920);
    tempMask1 = imdilate(tempMask1,ones(9));
    CorrMask1{corrCounter} = tempMask1;
    tempMask2 = poly2mask([Skeleton_x(:,2,corrCounter); Skeleton_x(1,2,corrCounter)],...
        [Skeleton_y(:,2,corrCounter);Skeleton_y(1,2,corrCounter)],1080,1920);
    tempMask2 = imdilate(tempMask2,ones(9));
    CorrMask2{corrCounter} = tempMask2;
end

% Update the fields in WS
WS.CorrMask1 = CorrMask1;
WS.CorrMask2 = CorrMask2;

end
