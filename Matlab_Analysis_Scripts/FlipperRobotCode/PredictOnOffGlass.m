% PREDICTONOFFGLASS Runs all crossing events through NN to predict glass/proper crossing
%
%  Uses the NN trained by Joseph and Baohua (netOnOffAmbig) to predict
%  whether a crossing event was a glass crossing (On), proper crossing (Off),
%  or ambiguous (Ambig). The output of the neural net is a probability
%  assigned to each of the three possibilities (summing to 1). The reason
%  the option of ambiguous is given to the NN is because it was found to
%  increase its classification performance. See the ROC/AUC figures for
%  reference on the relatively increase in performance.
%
%  The outputs of this function are later used by ClassifyOnOff to convert
%  from three probabilities to just one label (glass or proper crossing).
%  However, it is useful to break up the two steps in order to allow users
%  to experiment with different classification thresholds (as was done to
%  generate the ROC). The default threshold is 0.5. Further info can be
%  found in the ClassifyOnOff function regarding this threshold.

function WS = PredictOnOffGlass(WS)

% Port in the relevant fields from WS
NumCorridors = WS.NumCorridors;
NumGaps = WS.NumGaps;
inputFileName = WS.inputFileName;
indPos = WS.indPos;
indPosFrameBuffer = WS.indPosFrameBuffer;
FBFS = WS.FlipBinnedFlyStruct;
GapCentsXOdd = WS.GapCentsXOdd;
GapCentsYOdd = WS.GapCentsYOdd;
GapCentsXEven = WS.GapCentsXEven;
GapCentsYEven = WS.GapCentsYEven;
netOnOffAmbig = WS.netOnOffAmbig;

% Mean pixel value to normalize each frame to have
avgPixVal = 0.9;

% If inputFileName was saved with an absolute path for some reason, fix it
locOfAbsPathInStr = strfind(inputFileName,'All_Raw_Videos\');
if ~isempty(locOfAbsPathInStr)
    inputFileName = inputFileName(locOfAbsPathInStr+15:end);
end

% Initialize the video reader object to read in clips
reader1 = VideoReader(['..\0_Raw_Videos\',inputFileName]);

% Establish the size of the clips
frameWidth = 100;
frameHeight = 75;

% Which frames relative to the center frame to grab for the input
relativeFramesGrabbed = [30, 22, 18, 15, 12, 9, 3, 0, -3, -9, -12, -15, -18, -22, -30];
% Number of frames being fed in per sample
numFramesGrabbed = length(relativeFramesGrabbed);

% Number of flies being analyzed
numFlies = length(FBFS.ExpNum);

% Initialize various counters
i_th_flip = 1;
inOddFlipFrameCounter = 0;
inEvenFlipFrameCounter = 0;

% Initialize the matrix that holds all frames of a flip
numFramesInFlip = WS.flipRate*reader1.frameRate;
tempFrameHolderMat = zeros(reader1.height, reader1.width, numFramesInFlip);

% Now loop through the whole video
for globalFrameCounter = 1:reader1.NumFrames
    
    % The beginning and end frame numbers for each odd and even flip, but
    % check to make sure that these indices make sense
    if length(indPos) >= (4*(i_th_flip-1)+2)
        OddBeg_i    = indPos(4*(i_th_flip-1)+1)+indPosFrameBuffer;
        OddEnd_i    = indPos(4*(i_th_flip-1)+2)-indPosFrameBuffer;
    % If we're past the number of flips, just move forward to end of loop
    else
        continue
    end
    if length(indPos) >= (4*(i_th_flip-1)+4)
        EvenBeg_i   = indPos(4*(i_th_flip-1)+3)+indPosFrameBuffer;
        EvenEnd_i   = indPos(4*(i_th_flip-1)+4)-indPosFrameBuffer;
    % If we're past the number of flips, just move forward to end of loop
    else
        continue
    end
    
    % Read in the next frame
    currFrame = readFrame(reader1);
    
    % If the frame being read in belongs to an odd flip
    if (globalFrameCounter >= OddBeg_i) && (globalFrameCounter <= OddEnd_i)
        % Advance the oddFlipFrameCounter
        inOddFlipFrameCounter = inOddFlipFrameCounter + 1;
        % Grab all these frames and fill into matrix holding frames from i_th odd flip
        tempFrameHolderMat(:,:,inOddFlipFrameCounter) = currFrame(:,:,1);
    % If the frame being read in belongs to an even flip
    elseif (globalFrameCounter >= EvenBeg_i) && (globalFrameCounter <= EvenEnd_i)
        % Advance the evenFlipFrameCounter
        inEvenFlipFrameCounter = inEvenFlipFrameCounter + 1;
        % Grab all these frames and fill into matrix holding frames from i_th even flip
        tempFrameHolderMat(:,:,inEvenFlipFrameCounter) = currFrame(:,:,1);
    % If the frame being read in is the first frame after the last frame of an odd flip
    elseif inOddFlipFrameCounter
        % Determine the frame number of the first frame of the matrix
        firstFrameNumOfHolderMat = (globalFrameCounter - inOddFlipFrameCounter);
        % Reset the counter so that we don't get stuck in this loop later
        inOddFlipFrameCounter = 0;
        % Loop through all flies
        for flyCounter = 1:numFlies
            % Check that there are sufficiently many odd flips
            if length(FBFS.ExpNum(flyCounter).BehavData.OddFlips) >= i_th_flip
                % Loop through each gap size
                for gapCounter = 1:NumGaps
                    % Check that there are any valid up "crossings"
                    if FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesAbsFrameStart(gapCounter).GapID
                        % Grab the start and end frames for the "crossing(s)"
                        startFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesAbsFrameStart(gapCounter).GapID;
                        endFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesAbsFrameEnd(gapCounter).GapID;
                       
                        % Loop through all the "crossings" at that width in that flip for that fly
                        for eventCounter = 1:length(startFrames)
                            startFrame  = startFrames(eventCounter);
                            endFrame    = endFrames(eventCounter);

                            % Initialize the matrix that holds the time series of the crossing event
                            LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                            % Tracker for how many frames were written within sample
                            numFramesWritten = 0;
                            % Read in the appropriate frames of interest 
                            % (center frame + relativeFrames)
                            for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                                numFramesWritten = numFramesWritten + 1;
                                % Check to make sure that the frames requested are within the time points of the cross
                                % If not, fill in a white frame
                                if frameCounter < startFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                elseif frameCounter > endFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                % If frame requested is within time points of the cross
                                else
                                    % Read in the frame requested
                                    frame = tempFrameHolderMat(:,:,frameCounter-firstFrameNumOfHolderMat+1);
                                    % Only read in the pixels in the region of interest
                                    frame = frame((GapCentsYOdd(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYOdd(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                                 (GapCentsXOdd(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXOdd(gapCounter,flyCounter)+round(frameWidth/2)));
                                    % Normalize the frame
                                    normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                                end

                                % Fill up the matrix appropriately with the sample
                                % Invert the color for NN performance (255-normFrame)
                                LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                                    255-normFrame;
                                % Make sure to convert input to uint8 since NN is expecting that
                                LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                            end
                            
                            % Now run the NN on the crossing event
                            PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);

                            % Fill in the field with the NN's prediction probabilities
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(1);
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(2);
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(3);

                            % Print the flip counter to command window to give a measure of progress
%                             i_th_flip
                        end

                    % If there was no up "crossing" event, fill the probability fields with zeroes
                    else
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesAmbigProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesOffProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).UpCrossesOnProb(gapCounter).GapID = ...
                            0;
%                         i_th_flip
                    end
                    
                    % Check that there are any valid down "crossings"
                    if FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesAbsFrameStart(gapCounter).GapID
                        % Grab the start and end frames for the "crossing(s)"
                        startFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesAbsFrameStart(gapCounter).GapID;
                        endFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesAbsFrameEnd(gapCounter).GapID;
                        
                        % Loop through all the "crossings" at that width in that flip for that fly
                        for eventCounter = 1:length(startFrames)
                            startFrame  = startFrames(eventCounter);
                            endFrame    = endFrames(eventCounter);

                            % Initialize the matrix that holds the time series of the crossing event
                            LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                            % Tracker for how many frames were written within sample
                            numFramesWritten = 0;
                            % Read in the appropriate frames of interest 
                            % (center frame + relativeFrames)
                            for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                                numFramesWritten = numFramesWritten + 1;
                                % Check to make sure that the frames requested are within the time points of the cross
                                % If not, fill in a white frame
                                if frameCounter < startFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                elseif frameCounter > endFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                % If frame requested is within time points of the cross
                                else
                                    % Read in the frame requested
                                    frame = tempFrameHolderMat(:,:,frameCounter-firstFrameNumOfHolderMat+1);
                                    % Only read in the pixels in the region of interest
                                    frame = frame((GapCentsYOdd(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYOdd(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                                 (GapCentsXOdd(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXOdd(gapCounter,flyCounter)+round(frameWidth/2)));
                                    % Normalize the frame
                                    normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                                end

                                % Fill up the matrix appropriately with the sample
                                % Invert the color for NN performance (255-normFrame)
                                LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                                    255-normFrame;
                                % Make sure to convert input to uint8 since NN is expecting that
                                LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                            end
                                
                            % Now run the NN on the crossing event
                            PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);

                            % Fill in the field with the NN's prediction probabilities
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(1);
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(2);
                            FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(3);

                            % Print the flip counter to command window to give a measure of progress
%                             i_th_flip
                        end

                    % If there was no down "crossing" event, fill the probability fields with zeroes
                    else
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesAmbigProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesOffProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.OddFlips(i_th_flip).DownCrossesOnProb(gapCounter).GapID = ...
                            0;
%                         i_th_flip
                    end
                end
            end
        end
    % If the frame being read in is the first frame after the last frame of an odd flip
    elseif inEvenFlipFrameCounter
        % Determine the frame number of the first frame of the matrix
        firstFrameNumOfHolderMat = (globalFrameCounter - inEvenFlipFrameCounter);
        % Reset the counter so that we don't get stuck in this loop later
        inEvenFlipFrameCounter = 0;
        % Loop through all flies
        for flyCounter = 1:numFlies
            % Check that there are sufficiently many even flips
            if length(FBFS.ExpNum(flyCounter).BehavData.EvenFlips) >= i_th_flip
                % Loop through each gap size
                for gapCounter = 1:NumGaps
                    % Check that there are any valid up "crossings"
                    if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesAbsFrameStart(gapCounter).GapID
                        % Grab the start and end frames for the "crossing(s)"
                        startFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesAbsFrameStart(gapCounter).GapID;
                        endFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesAbsFrameEnd(gapCounter).GapID;
                        
                        % Loop through all the "crossings" at that width in that flip for that fly
                        for eventCounter = 1:length(startFrames)
                            startFrame  = startFrames(eventCounter);
                            endFrame    = endFrames(eventCounter);

                            % Initialize the matrix that holds the time series of the crossing event
                            LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                            % Tracker for how many frames were written within sample
                            numFramesWritten = 0;
                            % Read in the appropriate frames of interest 
                            % (center frame + relativeFrames)
                            for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                                numFramesWritten = numFramesWritten + 1;
                                % Check to make sure that the frames requested are within the time points of the cross
                                % If not, fill in a white frame
                                if frameCounter < startFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                elseif frameCounter > endFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                % If frame requested is within time points of the cross
                                else
                                    % Read in the frame requested
                                    frame = tempFrameHolderMat(:,:,frameCounter-firstFrameNumOfHolderMat+1);
                                    % Only read in the pixels in the region of interest
                                    frame = frame((GapCentsYEven(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYEven(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                                 (GapCentsXEven(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXEven(gapCounter,flyCounter)+round(frameWidth/2)));
                                    % Normalize the frame
                                    normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                                end

                                % Fill up the matrix appropriately with the sample
                                % Invert the color for NN performance (255-normFrame)
                                LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                                    255-normFrame;
                                % Make sure to convert input to uint8 since NN is expecting that
                                LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                            end
                            
                            % Now run the NN on the crossing event
                            PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);

                            % Fill in the field with the NN's prediction probabilities
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(1);
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(2);
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(3);

                            % Print the flip counter to command window to give a measure of progress
%                             i_th_flip
                        end

                    % If there was no up "crossing" event, fill the probability fields with zeroes
                    else
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesAmbigProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesOffProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).UpCrossesOnProb(gapCounter).GapID = ...
                            0;
%                         i_th_flip
                    end
                    
                    % Check that there are any valid down "crossings"
                    if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesAbsFrameStart(gapCounter).GapID
                        % Grab the start and end frames for the "crossing(s)"
                        startFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesAbsFrameStart(gapCounter).GapID;
                        endFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesAbsFrameEnd(gapCounter).GapID;
                        
                        % Loop through all the "crossings" at that width in that flip for that fly
                        for eventCounter = 1:length(startFrames)
                            startFrame  = startFrames(eventCounter);
                            endFrame    = endFrames(eventCounter);

                            % Initialize the matrix that holds the time series of the crossing event
                            LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                            % Tracker for how many frames were written within sample
                            numFramesWritten = 0;
                            % Read in the appropriate frames of interest 
                            % (center frame + relativeFrames)
                            for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                                numFramesWritten = numFramesWritten + 1;
                                % Check to make sure that the frames requested are within the time points of the cross
                                % If not, fill in a white frame
                                if frameCounter < startFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                elseif frameCounter > endFrame
                                    frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                                    normFrame = frame;
                                % If frame requested is within time points of the cross
                                else
                                    % Read in the frame requested
                                    frame = tempFrameHolderMat(:,:,frameCounter-firstFrameNumOfHolderMat+1);
                                    % Only read in the pixels in the region of interest
                                    frame = frame((GapCentsYEven(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYEven(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                                 (GapCentsXEven(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXEven(gapCounter,flyCounter)+round(frameWidth/2)));
                                    % Normalize the frame
                                    normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                                end

                                % Fill up the matrix appropriately with the sample
                                % Invert the color for NN performance (255-normFrame)
                                LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                                    255-normFrame;
                                % Make sure to convert input to uint8 since NN is expecting that
                                LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                            end
                                
                            % Now run the NN on the crossing event
                            PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);

                            % Fill in the field with the NN's prediction probabilities
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(1);
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(2);
                            FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                                PredAmbigOffOn(3);

                            % Print the flip counter to command window to give a measure of progress
%                             i_th_flip
                        end

                    % If there was no down "crossing" event, fill the probability fields with zeroes
                    else
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesAmbigProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesOffProb(gapCounter).GapID = ...
                            0;
                        FBFS.ExpNum(flyCounter).BehavData.EvenFlips(i_th_flip).DownCrossesOnProb(gapCounter).GapID = ...
                            0;
%                         i_th_flip
                    end
                end
            end
        end
        
        % Advance the flip ticker    
        i_th_flip = i_th_flip + 1;
        if mod(i_th_flip,10) == 0
            disp(['Data fed to NN up to flip ', num2str(i_th_flip), ' for file: ', inputFileName])
        end
            
    % The cases where the frame is not part of a flip nor the frame immediately following a flip
    else
        % Don't need to do anything except reset the counters
        inOddFlipFrameCounter = 0;
        inEvenFlipFrameCounter = 0;
    end
    
end

% Update the fields in WS
WS.FlipBinnedFlyStruct = FBFS;
WS.inputFileName = inputFileName;

end
