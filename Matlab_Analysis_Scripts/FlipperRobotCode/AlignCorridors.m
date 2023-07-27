% Takes all data in BehavData of a combined WS and produces a combined
% aligned WS in which all corridors have their coordinates transformed to
% be overlayed and perfectly straightened vertically. The new origin is at
% the center of the corridor, the y axis points up and x axis points right,
% and a set of x coordinates and orientations are provided in which the
% corridor is folded along its vertical axis of symmetry. Positive x in
% this coordinate system is defined as being into the gaps. All units of
% the aligned data are in real units (i.e., mm and sec).

% Meant to be used with the output of CombineExpStructs2

function WS = AlignCorridors(WS)

% Navigate to the data folder and establish where the appropriate save
% location is for the combined WS that was fed in
flipperRobotCodePath = pwd;
dataPath = [flipperRobotCodePath,'\..\..\Data\'];
cd(dataPath);
dataPath = pwd;
dataPath = [dataPath, '\'];
cd(flipperRobotCodePath);
% If the combined WS was composed of mismatching experiments, ask the user
% where the new aligned combined WS should be saved
if strcmpi(WS.directoryName, 'Potentially variable')
    GenotypeDirectory = uigetdir(dataPath,...
        'The experiments being aligned were combined from mismatching experiments. Select which folder you want to save the aligned WS to.');
% Otherwise, the aligned WS should be saved in the natural folder (which is
% listed as the WS's directoryName)
else
    GenotypeDirectory = [dataPath,WS.directoryName];
end
% Move to the folder in which we plan to save the aligned WS
cd(GenotypeDirectory);

% Conversions for pix to mm and frames to sec
PixelsPerMM = 15; % As measured from our videos
FramesPerSec = 30; % As recorded by our camera

% Go through and make the AlignedData field for each fly in WS
for flyCounter = 1:size(WS.FlipBinnedFlyStruct,2)
    WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips = ...
        WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips;
    WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips = ...
        WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips;
    % Figure out how many points compose the skeleton for a corridor (this
    % changes based on the number of gaps in a corridor, but it's typically
    % 36 because we typically use 4 gap cassettes)
    sizeSkel = size(WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x,1);
    % Compute the tilt angle (in degrees) of the corridors by computing the
    % angle between the vertical and the line that connects the bottom left
    % (skeleton point 1) and top left (skeleton point sizeSkel/2) points of
    % the corridor. We will use this angle later to rotate everything into
    % the vertical alignment.
    CorridorAngleOdd = atand((WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(1,1) - WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(sizeSkel/2,1))...
        / (WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(1,1) - WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(sizeSkel/2,1)));
    CorridorAngleEven = atand((WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(1,2) - WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(sizeSkel/2,2))...
        / (WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(1,2) - WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(sizeSkel/2,2)));
    % Compute the center of the corridor by taking the average of the top
    % left (skeleton point 1) and bottom right (skeleton point sizeSkel)
    % points of the corridor. We will use this center as our new origin.
    CorridorCenterXOdd  = 0.5*(WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(sizeSkel/2,1) + ...
                               WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(sizeSkel,1));
    CorridorCenterYOdd  = 0.5*(WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(sizeSkel/2,1) + ...
                               WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(sizeSkel,1));
    CorridorCenterXEven = 0.5*(WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(sizeSkel/2,2) + ...
                               WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(sizeSkel,2));
    CorridorCenterYEven = 0.5*(WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(sizeSkel/2,2) + ...
                               WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(sizeSkel,2));
    % Now we want to go through each even and odd flip and align all the data
    for oddFlipCounter = 1:size(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips,2)
        % If the flip had no data in it, do not transform the zeros
        if length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).CentroidX) == 1
            % If everything below is commented out up to the continue, then
            % flips with no data will contain no data (stored as [])
            % Otherwise, if the lines below are not commented out, then
            % flips with no data will contain single 0s in every field
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidX = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidY = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedVerticalTheta = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedCentroidX = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVerticalTheta = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelX = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelY = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelTheta = 0;
            continue
        % In a flip with data in it, transform into the aligned reference
        % frame by making the center of the corridor the origin and
        % rotating the entire corridor by the angle of the corridor and
        % then flipping the y axis to point upwards. We also convert all
        % the units into mm and sec (previously pixels and frames).
        else
            % Note that the x and y transformations are that of translating
            % to a new origin and then rotating COUNTERCLOCKWISE by 
            % CorridorAngle before then flipping the y axis
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).FrameInFlip/FramesPerSec;
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidX = ...
                ((WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).CentroidX - CorridorCenterXOdd)*cosd(CorridorAngleOdd) - ...
                 (WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).CentroidY - CorridorCenterYOdd)*sind(CorridorAngleOdd))/PixelsPerMM;
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidY = -1*...
                ((WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).CentroidY - CorridorCenterYOdd)*cosd(CorridorAngleOdd) + ...
                 (WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).CentroidX - CorridorCenterXOdd)*sind(CorridorAngleOdd))/PixelsPerMM;
            % To properly transform orientation of the ellipse, we rotate 
            % it by CorridorAngle and need to keep it bound from -90 to 90
            % degrees (which is the reason for the nested mods)
            % Also note the -1 corresponding to flipping the y axis in the
            % inner mod (which is why we add the corridor angle rather than
            % subtract it)
            % In this new AlignedOrientation, 0 degees corresponds to
            % moving perfectly in the direction of the corridor and +/-90
            % degrees corresponds to moving perfectly perpendicular to the
            % direction of the corridor (i.e., horizontal)
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedVerticalTheta = ...
                mod((mod(-1*WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).Orientation + CorridorAngleOdd + 90,180)-90),180)-90;
            % Now let's make the folded coordinates (i.e., the coordinates
            % in which the left half and right half of the corridors are
            % treated the same since there is a natural symmetry there)
            % This is easy to do in the new aligned coordinates because it
            % just consists of flipping negative x values to be positive
            % and changing the sign of the vertical theta whenever the x
            % values are negative.
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedCentroidX = ...
                abs(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidX);
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVerticalTheta = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedVerticalTheta.*...
                sign(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidX);
            % Finally, let's compute all the velocities (x, y, and theta)
            % in this new coordinate system
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelX = ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedCentroidX) ./ ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip);
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelY = ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidY) ./ ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip);
            % Doing this for theta is a little complicated because of the
            % periodicity of theta, but luckily we can use the unwrap
            % function in matlab to do this easily. Because we have an
            % unoriented ellipse, angles range from -90 to 90 instead of
            % -180 to 180 like the unwrap function expects, so to actually
            % get this function to work, we must multiply by 2 within the
            % function then divide by 2 outside of it. We also have to swap
            % back and forth from degrees and radians to make this work.
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelTheta = ...
                0.5*diff(rad2deg(unwrap(deg2rad(2*WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedVerticalTheta)))) ./ ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip);
        end
    end
    for evenFlipCounter = 1:size(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips,2)
        % If the flip had no data in it, do not transform the zeros
        if length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).CentroidX) == 1
            % If everything below is commented out up to the continue, then
            % flips with no data will contain no data (stored as [])
            % Otherwise, if the lines below are not commented out, then
            % flips with no data will contain single 0s in every field
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidX = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidY = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedVerticalTheta = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedCentroidX = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVerticalTheta = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelX = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelY = 0;
%             WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelTheta = 0;
            continue
        % In a flip with data in it, transform into the aligned reference
        % frame by making the center of the corridor the origin and
        % rotating the entire corridor by the angle of the corridor and
        % then flipping the y axis to point upwards. We also convert all
        % the units into mm and sec (previously pixels and frames).
        else
            % Note that the x and y transformations are that of translating
            % to a new origin and then rotating COUNTERCLOCKWISE by 
            % CorridorAngle before then flipping the y axis
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedTimeInFlip = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).FrameInFlip/FramesPerSec;
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidX = ...
                ((WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).CentroidX - CorridorCenterXEven)*cosd(CorridorAngleEven) - ...
                 (WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).CentroidY - CorridorCenterYEven)*sind(CorridorAngleEven))/PixelsPerMM;
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidY = -1*...
                ((WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).CentroidY - CorridorCenterYEven)*cosd(CorridorAngleEven) + ...
                 (WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).CentroidX - CorridorCenterXEven)*sind(CorridorAngleEven))/PixelsPerMM;
            % To properly transform orientation of the ellipse, we rotate 
            % it by CorridorAngle and need to keep it bound from -90 to 90
            % degrees (which is the reason for the nested mods)
            % Also note the -1 corresponding to flipping the y axis in the
            % inner mod (which is why we add the corridor angle rather than
            % subtract it)
            % In this new AlignedOrientation, 0 degees corresponds to
            % moving perfectly in the direction of the corridor and +/-90
            % degrees corresponds to moving perfectly perpendicular to the
            % direction of the corridor (i.e., horizontal)
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedVerticalTheta = ...
                mod((mod(-1*WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).Orientation + CorridorAngleEven + 90,180)-90),180)-90;
            % Now let's make the folded coordinates (i.e., the coordinates
            % in which the left half and right half of the corridors are
            % treated the same since there is a natural symmetry there)
            % This is easy to do in the new aligned coordinates because it
            % just consists of flipping negative x values to be positive
            % and changing the sign of the vertical theta whenever the x
            % values are negative.
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedCentroidX = ...
                abs(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidX);
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVerticalTheta = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedVerticalTheta.*...
                sign(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidX);
            % Finally, let's compute all the velocities (x, y, and theta)
            % in this new coordinate system
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelX = ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedCentroidX) ./ ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedTimeInFlip);
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelY = ...=
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidY) ./ ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedTimeInFlip);
            % Doing this for theta is a little complicated because of the
            % periodicity of theta, but luckily we can use the unwrap
            % function in matlab to do this easily. Because we have an
            % unoriented ellipse, angles range from -90 to 90 instead of
            % -180 to 180 like the unwrap function expects, so to actually
            % get this function to work, we must multiply by 2 within the
            % function then divide by 2 outside of it. We also have to swap
            % back and forth from degrees and radians to make this work.
            WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelTheta = ...=
                0.5*diff(rad2deg(unwrap(deg2rad(2*WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedVerticalTheta)))) ./ ...
                diff(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedTimeInFlip);
        end
    end
    % Now just make an aligned skeleton by performing the same transformations
    % AlignedSkeleton_x for odd flips
    WS.FlipBinnedFlyStruct(flyCounter).IdData.AlignedSkeleton_x(:,1) = ...
        ((WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(:,1) - CorridorCenterXOdd)*cosd(CorridorAngleOdd) - ...
         (WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(:,1) - CorridorCenterYOdd)*sind(CorridorAngleOdd))/PixelsPerMM;
    % AlignedSkeleton_y for odd flips
    WS.FlipBinnedFlyStruct(flyCounter).IdData.AlignedSkeleton_y(:,1) = -1*...
        ((WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(:,1) - CorridorCenterYOdd)*cosd(CorridorAngleOdd) + ...
         (WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(:,1) - CorridorCenterXOdd)*sind(CorridorAngleOdd))/PixelsPerMM;
    % AlignedSkeleton_x for even flips
    WS.FlipBinnedFlyStruct(flyCounter).IdData.AlignedSkeleton_x(:,2) = ...
        ((WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(:,2) - CorridorCenterXEven)*cosd(CorridorAngleEven) - ...
         (WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(:,2) - CorridorCenterYEven)*sind(CorridorAngleEven))/PixelsPerMM;
    % AlignedSkeleton_y for even flips
    WS.FlipBinnedFlyStruct(flyCounter).IdData.AlignedSkeleton_y(:,2) = -1*...
        ((WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_y(:,2) - CorridorCenterYEven)*cosd(CorridorAngleEven) + ...
         (WS.FlipBinnedFlyStruct(flyCounter).IdData.Skeleton_x(:,2) - CorridorCenterXEven)*sind(CorridorAngleEven))/PixelsPerMM;
end

% Now save WS with the alignment completed
% Ask the user where WS should  be saved if it came from a mismatched
% set of combined experiments
if strcmpi(WS.directoryName, 'Potentially variable')
    combinedWSName = input('What should this combined and aligned WS be named within your previously selected folder?\n','s');
    save([combinedWSName,'_WS_Combined.mat'],'WS');
% If the WS has a natural directory, just save it there as [...]WS_Combined
else
    save([WS.directoryName,'_WS_Combined.mat'],'WS')
end

% Navigate back to the directory with all the code
cd(flipperRobotCodePath);

end