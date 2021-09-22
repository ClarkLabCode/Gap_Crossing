NumCorridors = 7;
NumGaps = 4;

% Mean pixel value to normalize each frame to have
avgPixVal = 0.9;

% Initialize the video reader object to read in clips
reader = VideoReader('2021-07-07 14-11-22_IsoD1_30C_berberine_chloride_40_min_training.mp4');

% Establish the size of the clips
frameWidth = 100;
frameHeight = 75;

% Which frames relative to the center frame to grab for the input
relativeFramesGrabbed = [30, 22, 18, 15, 12, 9, 3, 0, -3, -9, -12, -15, -18, -22, -30];
% Number of frames being fed in per sample
numFramesGrabbed = length(relativeFramesGrabbed);

% Matrix that holds each sample. The direction of time is upwards in the image
LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1),length(find(CrossFrameIDLabels)));
% Vector that holds the labels for each sample
LabeledCrossingVecLabel = zeros(length(find(CrossFrameIDLabels)),1);

% Initialize the counter for the sample number
sampleCounter = 0;

% Go through the loop for each fly at each width for each orientation
% (even/odd flip) for each event number
for eventCounter = 1:length(CrossFrameIDLabels(1,1,:))
    for gapCounter = 1:NumGaps
        for flyCounter = 1:NumCorridors
            % Check that there is indeed this many events for a given fly
            % at a given width for odd flips
            if CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) ~= 0
                % No semi-colon so that progress can be tracked
                sampleCounter = sampleCounter + 1
                % Tracker for how many frames were written within sample
                numFramesWritten = 0;
                % Read in the appropriate frames of interest 
                % (center frame + relativeFrames)
                for frameCounter = ...
                (round(((CrossFrameIDMat(flyCounter,gapCounter,1,eventCounter))+(CrossFrameIDMat(flyCounter,gapCounter,2,eventCounter)))/2) + ...
                relativeFramesGrabbed)
                    numFramesWritten = numFramesWritten + 1;
                    % Check to make sure that the frames requested are within the time points of the flip
                    % If not, fill in a white frame
                    if frameCounter < (CrossFrameIDMat(flyCounter,gapCounter,1,eventCounter))
                        frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                        normFrame = frame;
                    elseif frameCounter > (CrossFrameIDMat(flyCounter,gapCounter,2,eventCounter))
                        frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                        normFrame = frame;
                    else
                        frame = read(reader, frameCounter);
                        % Only reads in one channel since not RGB
                        frame = frame(:,:,1);
                        % Only read in the pixels in the region of interest
                        frame = frame((y_cents_1(gapCounter,flyCounter)-round(frameHeight/2)):(y_cents_1(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                     (x_cents_1(gapCounter,flyCounter)-round(frameWidth/2)):(x_cents_1(gapCounter,flyCounter)+round(frameWidth/2)));
                        % Normalize the frame
                        normFrame = double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                    end
                    % Fill up the matrix appropriately with the sample
                    LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:,sampleCounter) = ...
                        1-normFrame;
                    % Fill up the vector appropriately with the sample label
                    LabeledCrossingVecLabel(sampleCounter) = CrossFrameIDLabels(flyCounter,gapCounter,eventCounter);
                end
            end
                                
                
            
            % Check that there is indeed this many events for a given fly
            % at a given width for even flips
            if CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) ~= 0
                sampleCounter = sampleCounter + 1;
                % Tracker for how many frames were written within sample
                numFramesWritten = 0;
                % Read in the appropriate frames of interest 
                % (center frame + relativeFrames)
                for frameCounter = ...
                (round(((CrossFrameIDMat(flyCounter,gapCounter+NumGaps,1,eventCounter))+(CrossFrameIDMat(flyCounter,gapCounter+NumGaps,2,eventCounter)))/2) + ...
                relativeFramesGrabbed)
                    numFramesWritten = numFramesWritten + 1;
                    % Check to make sure that the frames requested are within the time points of the flip
                    % If not, fill in a white frame
                    if frameCounter < CrossFrameIDMat(flyCounter,gapCounter+NumGaps,1,eventCounter)
                        frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                        normFrame = frame;
                    elseif frameCounter > CrossFrameIDMat(flyCounter,gapCounter+NumGaps,2,eventCounter)
                        frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                        normFrame = frame;
                    else
                        frame = read(reader, frameCounter);
                        % Only reads in one channel since not RGB
                        frame = frame(:,:,1);
                        % Only read in the pixels in the region of interest
                        frame = frame((y_cents_2(gapCounter,flyCounter)-round(frameHeight/2)):(y_cents_2(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                     (x_cents_2(gapCounter,flyCounter)-round(frameWidth/2)):(x_cents_2(gapCounter,flyCounter)+round(frameWidth/2)));
                        % Normalize the frame
                        normFrame = double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                    end
                    % Fill up the matrix appropriately with the sample
                    LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:,sampleCounter) = ...
                        1-normFrame;
                    LabeledCrossingVecLabel(sampleCounter) = CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter);
                end
            end
        end
    end
end

% % Rescales to grayscale
% LabeledCrossingMatInput = LabeledCrossingMatInput/255;



% LabeledCrossingMatInputAllFlips = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1),4*length(find(CrossFrameIDLabels)));
% 
% LabeledCrossingMatInputAllFlips(:,:,1:end/4) = LabeledCrossingMatInput;
% LabeledCrossingMatInputAllFlips(:,:,(end/4+1):end/2) = LabeledCrossingMatInput2;
% LabeledCrossingMatInputAllFlips(:,:,(end/2+1):end*3/4) = LabeledCrossingMatInput3;
% LabeledCrossingMatInputAllFlips(:,:,(end*3/4+1):end) = LabeledCrossingMatInput4;
% 
% LabeledCrossingVecLabelAllFlips = repmat(LabeledCrossingVecLabel,[4 1]);


% % Randomly shuffle the samples
% shuffleOrder = randperm(length(LabeledCrossingMatInput(1,1,:)));
% shuffledInputMat = LabeledCrossingMatInput(:,:,shuffleOrder);
% shuffledLabelVec = LabeledCrossingVecLabel(shuffleOrder);

% % Normalize the grayscale per sample to go from 0 to 1
% LabeledCrossingMatInputRenorm = LabeledCrossingMatInput./(max(max(LabeledCrossingMatInput,[],1),[],2));
% LabeledCrossingMatInput = LabeledCrossingMatInputRenorm;

% % Get the dimensions of LabeledCrossingMatInput
% [y, x, numSamples] = size(LabeledCrossingMatInput);
% 
% % Normalize the grayscale per sample based on the average pixel value
% LabeledCrossingMatInputRenorm = LabeledCrossingMatInput./(sum(sum(LabeledCrossingMatInput,2),1))*y*x*avgPixVal;
% LabeledCrossingMatInput = LabeledCrossingMatInputRenorm;