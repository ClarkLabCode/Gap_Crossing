% Figures out which frames in the video correspond to flips of cassette

% Inputs:
% inputFileName     = Name of video file to be analyzed, must include extension
% flipRate          = Frequency of flipping (in seconds)
% plotYN            = Whether or not to plot the change vectors (1 for yes)

% Outputs:
% indPos            = Index of frame at which flips happen
% reader1           = VideoReader object that reads the inputFileName
% frameMarkerVec    = Vector that is 1 only where subsequent frames are sufficiently different
% changeVecDil      = Vector used to find frameMarkerVec, see code for more info
% dotProdVec        = Vector of dot products of subsequent frames, used to find
%                     the above vectors

function [indPos, reader1, frameMarkerVec, changeVecDil, dotProdVec] = ...
    FindFramesForFlip(inputFileName, flipRate, plotYN)

reader1 = VideoReader(inputFileName);
frameCount = floor(reader1.Duration * reader1.FrameRate) - 1;
resolution = [reader1.Width reader1.Height];
dotProdVec = zeros(frameCount,1);
changeVec = zeros(frameCount-3,1); % 3 frame rolling tolerance
frameMarkerVec = zeros(frameCount-4,1);
numFramePerFlip = (flipRate-2)*reader1.FrameRate;   % giving 2 sec tolerance
vid_array = zeros(resolution(2), resolution(1), numFramePerFlip, 'uint8');

% Grab all frames within the first flip and put them into vid_array
for i = 1:numFramePerFlip
    vid_array(:,:,i) = rgb2gray(readFrame(reader1));
end

% Find the bg frame within first flip (dimension 3 is time in vid_array)
bg_img1 = uint8(mean(vid_array,3));

% Figure out when subsequent frames are very different by computing the
% dot product of each frame with the bg frame and integrating in both dims
for i = 1:numFramePerFlip
    % Compute the dot product of the bg img with each frame and integrate
    dotProdVec(i) = sum(sum(bg_img1.*vid_array(:,:,i)));
    % Used to skip the first 3 frames for indexing reasons
    if i<=3
        continue
    end
    % Creates changeVec which is 1 if a frame is within +/- 5% of the
    % integrated dot product of its previous 3 frames and 0 otherwise
    changeVec(i-3) = (dotProdVec(i)>0.95*dotProdVec(i-1))*(dotProdVec(i)<1.05*dotProdVec(i-1))...
    *(dotProdVec(i)>0.95*dotProdVec(i-2))*(dotProdVec(i)<1.05*dotProdVec(i-2))...
    *(dotProdVec(i)>0.95*dotProdVec(i-3))*(dotProdVec(i)<1.05*dotProdVec(i-3));
end

% Repeats above loop but for after the first flip
for i = (numFramePerFlip+1):frameCount
    dotProdVec(i) = sum(sum(bg_img1.*rgb2gray(readFrame(reader1))));
    changeVec(i-3) = (dotProdVec(i)>0.95*dotProdVec(i-1))*(dotProdVec(i)<1.05*dotProdVec(i-1))...
    *(dotProdVec(i)>0.95*dotProdVec(i-2))*(dotProdVec(i)<1.05*dotProdVec(i-2))...
    *(dotProdVec(i)>0.95*dotProdVec(i-3))*(dotProdVec(i)<1.05*dotProdVec(i-3));
end

% Using circshift to dilate out the 0s within changeVec 
% e.g. if changeVec = [0 1 1 1 0 1 0 1 1 1]
% then changeVecDil = [0 0 1 0 0 0 0 0 1 1]
% Do this to err on the side of caution for flip marking purposes
changeVecFor = circshift(changeVec,1);
changeVecBack = circshift(changeVec,-1);
changeVecDil = changeVec.*changeVecFor.*changeVecBack;

% Create frameMarkerVec which is 1 whenever subsequent frames differ and
% 0 otherwise
for i = 1:(frameCount-3)
    if i == 1
        continue
    end
    frameMarkerVec(i-1) = ((changeVecDil(i-1)-changeVecDil(i)) ~= 0);
end

% Using frameMarkerVec, find the corresponding frames in which a flip of
% the robot happens and track these frame numbers in indPos
indPos = find(frameMarkerVec);
indPos = indPos + 4; % Unshifting the frame count

% Currently indPos doesn't have the first frame of the video and may not
% have the last frame of the video. Check if it has the last frame of the
% video by seeing whether it has an even or odd number of elements. If it
% has the last frames of a video, then it will have an odd number, so no
% need to add the last frame of the video as the last element
if ((-1)^(length(indPos)) == 1)
    indPos(end+1) = frameCount; % putting end points on start/end frame
end        

% Add the first frame of the video as the first element of indPos
indPos(end+1) = 1;
indPos = circshift(indPos,1);

% If plotting is requested
if plotYN == 1

    % Plot changeVecDil and dotProdVec just to see what they look like
    % If it's working as expected, then changeVecDil should be 0 the vast
    % majority of the time and 1 periodically when the flip happens
    % dotProdVec should look like two different baseline levels separated
    % by sudden changes during flip transitions
    figure(1)
    plot(changeVecDil);
    figure(2)
    plot(dotProdVec);

end

end