% Extracts all the desired raw info from the video and puts into struct

% Inputs:
% inputFileName     = Name of video file to be analyzed, must include extension
% directoryName     = Base name of the file where info will be saved
% sizeThreshCutOff  = Min number of pixels to count as fly, ~100 is good
% indPosFrameBuffer = How many frames before and after flip to skip, ~5 is good
% erodePix          = How many times to run "imerode" on bg sub frames, usually 1

% Outputs:
% finalStats    = Structure that holds all info from the video, each row
%                 corresponds to a particular fly at a particular time
% AreaVec       = Vector of areas used as intermediate for removing non-flies
% AreaVecLog    = Logical vector used to threshold AreaVec

function [finalStats, AreaVec, AreaVecLog] = ...
    RawFlyInfoExtracter(inputFileName, directoryName, sizeThreshCutOff, indPosFrameBuffer, erodePix)

% Open a video reader and grab the important quantities from it
reader1 = VideoReader(inputFileName);
frameCount = floor(reader1.Duration * reader1.FrameRate) - 1;
resolution = [reader1.Width reader1.Height];

numFlips = (length(indPos)/2); % How many flips take place in a video

% Go through frames of one flip at a time and perform bg subtraction, then
% threshold the frames (keeping pixels that are 1/5th the intensity of the
% max pixel), then erode the image (helps remove streaks from cassettes),
% then run regionprops on each frame and save all the relevant info
for flipCounter = 1:numFlips
    % tempArray holds all frames of the flip (3rd dim is RGB, so no use)
    tempArray = read(reader1, [indPos(2*flipCounter-1)+indPosFrameBuffer indPos(2*flipCounter)-indPosFrameBuffer]);
    sizeOfTemp = size(tempArray);
    flipArray = zeros(sizeOfTemp(1), sizeOfTemp(2), sizeOfTemp(4));
    
    % flipArray converts tempArray to grayscale by removing 3rd dim
    flipArray(:,:,:) = tempArray(:,:,1,:);
    
    % Save the bg image for each flip so we can use later if needed
    bg_img_name = strcat(directoryName,'_bg_flip',num2str(flipCounter),'.png');
    imwrite(uint8(mean(flipArray,3)), bg_img_name);
    
    % Bg subtract flipArray, then threshold
    flipArray = mean(flipArray,3) - flipArray;
    flipArrayBW = flipArray > (max(max(flipArray))/5);
    
    % Erode flipArray by using the circshift trick below
    flipArrayBWErode = flipArrayBW;
    for LRShift = -erodePix:erodePix
        for UDShift = -erodePix:erodePix
            flipArrayBWErode = flipArrayBWErode.*circshift(flipArrayBW, [LRShift UDShift]);
        end
    end
    flipArrayBWErode = logical(flipArrayBWErode);
    
    % Need to separate out first flip from rest of flips so that we can
    % create the structure finalStats appropriately the first time
    % Each row of finalStats corresponds to a fly in a particular frame and
    % exists for all flies across all frames of the video. This structure
    % will later be converted to a fly-centric structure where each row
    % corresponds solely to a fly and contains all temporal info within it
    if flipCounter == 1
        for frameCountInFlip = 1:sizeOfTemp(4)
            stats = regionprops('struct', flipArrayBWErode(:,:,frameCountInFlip),...
                'Area', 'Centroid', 'BoundingBox', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
            for flyCounter = 1:length(stats)
                stats(flyCounter).FlipNumber = flipCounter;
                stats(flyCounter).FrameInFlip = frameCountInFlip;
                stats(flyCounter).AbsoluteTime = indPos(2*flipCounter-1)+indPosFrameBuffer + frameCountInFlip - 1;
            end
            if frameCountInFlip == 1
                finalStats = stats;
            else
                finalStats(end+1:end+length(stats)) = stats;
            end
        end
    else
        for frameCountInFlip = 1:sizeOfTemp(4)
            stats = regionprops('struct', flipArrayBWErode(:,:,frameCountInFlip),...
                'Area', 'Centroid', 'BoundingBox', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
            for flyCounter = 1:length(stats)
                stats(flyCounter).FlipNumber = flipCounter;
                stats(flyCounter).FrameInFlip = frameCountInFlip;
                stats(flyCounter).AbsoluteTime = indPos(2*flipCounter-1)+indPosFrameBuffer + frameCountInFlip - 1;
            end
            finalStats(end+1:end+length(stats)) = stats;
        end
    end
end

% Remove all "flies" which do not meet the min fly size threshold
AreaVec = [finalStats(:).Area];
AreaVecLog = AreaVec < sizeThreshCutOff;
finalStats(AreaVecLog) = [];

end