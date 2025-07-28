% SKELETONTOCORRMASKS Takes user-made skeleton to generate corridor masks
% 
%  Once the user has made the skeleton (CorridorSkeletonFinder), this
%  function takes the labeled skeleton points and computes the mask for
%  each individual corridor. For a reminder of how a corridor is defined,
%  please refer to the diagram in the README.

function WS = SkeletonToCorrMasks(WS)

% Port in the relevant fields from WS
NumGaps = WS.NumGaps;
NumCorridors = WS.NumCorridors;
Skeleton_x = WS.Skeleton_x;
Skeleton_y = WS.Skeleton_y;
GapOrientation = WS.GapOrientation;
% Grab the resolution of the frames so we can make appropriately sized masks
if isfield(WS,'resolution')
    resolution = WS.resolution;
% To make this backwards compatible with data that was processed before
% this resolution was saved, we add this clause that auto-fills the
% resolution if the gap orientation was saved but not the resolution
else
    if strcmpi(GapOrientation,'Vertical')
        resolution = [1920,1080];
    elseif strcmpi(GapOrientation,'Horizontal')
        resolution = [1080,1920];
    else
        error('GapOrientation entered was unexpected. Please use "Vertical" or "Horizontal".')
    end
end


% Initialize the mask cells
CorrMask1 = cell(1, NumCorridors);
CorrMask2 = cell(1, NumCorridors);

% Initialize a temporary mask to be used
% Make sure to do this in a way that allows for vertical or horizontal gaps
tempMask1 = false(resolution(2),resolution(1));
tempMask2 = false(resolution(2),resolution(1));

% Loop through the skeleton data for each corridor and fill in for each orientation
% For poly2mask to work, we need to close the polygon by repeating the first point
% We then dilate the mask slightly to make sure we capture the full corridor
for corrCounter = 1:NumCorridors
    tempMask1 = poly2mask([Skeleton_x(:,1,corrCounter); Skeleton_x(1,1,corrCounter)],...
        [Skeleton_y(:,1,corrCounter);Skeleton_y(1,1,corrCounter)],resolution(2),resolution(1));
    tempMask1 = imdilate(tempMask1,ones(9));
    CorrMask1{corrCounter} = tempMask1;
    tempMask2 = poly2mask([Skeleton_x(:,2,corrCounter); Skeleton_x(1,2,corrCounter)],...
        [Skeleton_y(:,2,corrCounter);Skeleton_y(1,2,corrCounter)],resolution(2),resolution(1));
    tempMask2 = imdilate(tempMask2,ones(9));
    CorrMask2{corrCounter} = tempMask2;
end

% Update the fields in WS
WS.CorrMask1 = CorrMask1;
WS.CorrMask2 = CorrMask2;

end
