% FINDFRAMESFORFLIP Figures out which video frames correspond to flips of cassette
%
%  This is done by looking for frame sequences in which there is substantial
%  frame-to-frame differences and for multiple consecutive frames. This is a 
%  signature of a flip of the cassette because as the cassette flips, it no
%  longer blocks as much of the backlighting IR LEDs, so there is an overall
%  large change in pixel-to-pixel luminance.

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

function WS = FindFramesForFlip(WS)

% function [indPos, reader1, frameMarkerVec, changeVecDil, dotProdVec] = ...
%     FindFramesForFlip(inputFileName, flipRate, plotYN)

% Port in the relevant fields from WS
inputFileName       = WS.inputFileName;
flipRate            = WS.flipRate;
indPosFrameBuffer   = WS.indPosFrameBuffer;
plotYN              = 0;

reader1 = VideoReader(['..\0_Raw_Videos\',inputFileName]);
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
    changeVec(i-3) = (dotProdVec(i)>0.99*dotProdVec(i-1))*(dotProdVec(i)<1.01*dotProdVec(i-1))...
    *(dotProdVec(i)>0.99*dotProdVec(i-2))*(dotProdVec(i)<1.01*dotProdVec(i-2))...
    *(dotProdVec(i)>0.99*dotProdVec(i-3))*(dotProdVec(i)<1.01*dotProdVec(i-3));
end

% Repeats above loop but for after the first flip
for i = (numFramePerFlip+1):frameCount
    dotProdVec(i) = sum(sum(bg_img1.*rgb2gray(readFrame(reader1))));
    changeVec(i-3) = (dotProdVec(i)>0.99*dotProdVec(i-1))*(dotProdVec(i)<1.01*dotProdVec(i-1))...
    *(dotProdVec(i)>0.99*dotProdVec(i-2))*(dotProdVec(i)<1.01*dotProdVec(i-2))...
    *(dotProdVec(i)>0.99*dotProdVec(i-3))*(dotProdVec(i)<1.01*dotProdVec(i-3));
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

% See what the spacing (in frames) between subsequent marked frames is
diffOfIndPos = diff(indPos);

% Estimate the number of frames it takes for a flip to happen
% Do this by finding the smallest spacing between frames in the beginning
% of indPos (skipping the first entry) then taking the max of that number
% and the one that follows it by 2 entries (helps avoid errors this way)
% [~, locOfMin] = min(diffOfIndPos(2:5));
% numFramesForFlipTime = max(diffOfIndPos(locOfMin+1),diffOfIndPos(locOfMin+3));
numFramesForFlipTime = 25;

% Estimate the number of frames that are contained in a flip
numFramesWithinFlip = flipRate*reader1.FrameRate-3;

% Sometimes, the frames within one flip get counted as frames from several
% different flips, so this is a way to resolve that and group them all into
% one flip
% First look for which entries in indPos definitely belong to stationary periods
IndPosWithinFlip = diffOfIndPos > numFramesWithinFlip/(flipRate/2);
% Find where the array goes from non-zero to zero and vice versa
transitions = diff([0; IndPosWithinFlip == 0; 0]); 
% Mark these locations
runstarts = find(transitions == 1);
runends = find(transitions == -1);
runlengths = runends - runstarts;
% Look for regions in the vector where there are more than 1 consecutive
% indices that correspond to non-stationary periods and only keep those
runstarts(runlengths < 2) = [];
runends(runlengths < 2) = [];
% Create an array that holds the limits of these runs
runLims = [runstarts, runends];
% Initialize a vector that contains which indices to exclude from indPos
RunExcVec = [];
% Fill this array with the frames between each run start and end
for i = 1:size(runLims,1)
    RunExcVec = [RunExcVec, (runLims(i, 1)+1):(runLims(i,2)-1)];
end
% Now duplicate indPos into indPosPrime
indPosPrime = indPos;
% Finally, remove the entries in indPosPrime that are incorrect in indPos
indPosPrime(RunExcVec) = [];

% Check that indPosPrime is indeed fixed
diffOfIndPosPrime = diff(indPosPrime);
diffOfIndPosPrimeWithinExpectedRange = ...
(diffOfIndPosPrime < numFramesForFlipTime+15)&(diffOfIndPosPrime > numFramesForFlipTime-10)|...
(diffOfIndPosPrime < numFramesWithinFlip+10)&(diffOfIndPosPrime > numFramesWithinFlip-15);
% The first and last entries are always going to be 0 but aren't problems,
% so just set them to 1
diffOfIndPosPrimeWithinExpectedRange(1) = 1;
diffOfIndPosPrimeWithinExpectedRange(end) = 1;

% If indPosPrime isn't fixed, throw an error
if sum(~diffOfIndPosPrimeWithinExpectedRange)
    error('Something is wrong with finding the frames in which flips took place.')
end

% Replace indPos with indPosPrime since we now know it to be fixed
indPos = indPosPrime;

% Currently indPos may not have the first frame of the video and may not
% have the last frame of the video. Check if it has the first frame. If
% not, then add the first frame. Then check if it has the last frame of the
% video by seeing whether it has an even or odd number of elements. If it
% doesn't have the last frame, then it will have an odd number of elements, 
% so we need to add the last frame of the video as the last element

% Add the first frame of the video as the first element of indPos if it
% isn't already included
if indPos(1) > (1+indPosFrameBuffer)
    indPos(end+1) = 1;
    indPos = circshift(indPos,1);
end

% Check if the last frame is included, and if it isn't, then add to indPos
if ((-1)^(length(indPos)) ~= 1)
    indPos(end+1) = frameCount; % putting end points on start/end frame
end        

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

% Update the fields in WS
WS.indPos               = indPos;
WS.frameMarkerVec       = frameMarkerVec;
WS.changeVecDil         = changeVecDil;
WS.dotProdVec           = dotProdVec;

end
