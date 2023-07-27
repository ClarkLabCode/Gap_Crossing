function WS = PredictOnOffGlass(WS)

NumCorridors = WS.NumCorridors;
NumGaps = WS.NumGaps;
inputFileName = WS.inputFileName;
FBFS = WS.FlipBinnedFlyStruct;
GapCentsXOdd = WS.GapCentsXOdd;
GapCentsYOdd = WS.GapCentsYOdd;
GapCentsXEven = WS.GapCentsXEven;
GapCentsYEven = WS.GapCentsYEven;
netOnOffAmbig = WS.netOnOffAmbig;

% Mean pixel value to normalize each frame to have
avgPixVal = 0.9;

% Initialize the video reader object to read in clips
reader1 = VideoReader(inputFileName);

% Establish the size of the clips
frameWidth = 100;
frameHeight = 75;

% Which frames relative to the center frame to grab for the input
relativeFramesGrabbed = [30, 22, 18, 15, 12, 9, 3, 0, -3, -9, -12, -15, -18, -22, -30];
% Number of frames being fed in per sample
numFramesGrabbed = length(relativeFramesGrabbed);

% Go through the loop for each fly at each width for each orientation
% (even/odd flip) for each event number
for flyCounter = 1:length(FBFS.ExpNum)
    for oddFlipCounter = 1:length(FBFS.ExpNum(flyCounter).BehavData.OddFlips)
        for gapCounter = 1:NumGaps
            if FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesAbsFrameStart(gapCounter).GapID
                startFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesAbsFrameStart(gapCounter).GapID;
                endFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesAbsFrameEnd(gapCounter).GapID;
                
                for eventCounter = 1:length(startFrames)
                    startFrame  = startFrames(eventCounter);
                    endFrame    = endFrames(eventCounter);
                    
                    LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                    % Tracker for how many frames were written within sample
                    numFramesWritten = 0;
                    % Read in the appropriate frames of interest 
                    % (center frame + relativeFrames)
                    for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                        numFramesWritten = numFramesWritten + 1;
                        % Check to make sure that the frames requested are within the time points of the flip
                        % If not, fill in a white frame
                        if frameCounter < startFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        elseif frameCounter > endFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        else
                            frame = read(reader1, frameCounter);
                            % Only reads in one channel since not RGB
                            frame = frame(:,:,1);
                            % Only read in the pixels in the region of interest
                            frame = frame((GapCentsYOdd(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYOdd(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                         (GapCentsXOdd(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXOdd(gapCounter,flyCounter)+round(frameWidth/2)));
                            % Normalize the frame
                            normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                        end
                        % Fill up the matrix appropriately with the sample
                        LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                            255-normFrame;
                        LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                    end
                    PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(1);
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(2);
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(3);
                    oddFlipCounter
                end
            else
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesAmbigProb(gapCounter).GapID = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesOffProb(gapCounter).GapID = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCrossesOnProb(gapCounter).GapID = ...
                    0;
                oddFlipCounter
            end
            
            if FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesAbsFrameStart(gapCounter).GapID
                startFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesAbsFrameStart(gapCounter).GapID;
                endFrames = FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesAbsFrameEnd(gapCounter).GapID;
                
                for eventCounter = 1:length(startFrames)
                    startFrame  = startFrames(eventCounter);
                    endFrame    = endFrames(eventCounter);
                    
                    LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                    % Tracker for how many frames were written within sample
                    numFramesWritten = 0;
                    % Read in the appropriate frames of interest 
                    % (center frame + relativeFrames)
                    for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                        numFramesWritten = numFramesWritten + 1;
                        % Check to make sure that the frames requested are within the time points of the flip
                        % If not, fill in a white frame
                        if frameCounter < startFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        elseif frameCounter > endFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        else
                            frame = read(reader1, frameCounter);
                            % Only reads in one channel since not RGB
                            frame = frame(:,:,1);
                            % Only read in the pixels in the region of interest
                            frame = frame((GapCentsYOdd(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYOdd(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                         (GapCentsXOdd(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXOdd(gapCounter,flyCounter)+round(frameWidth/2)));
                            % Normalize the frame
                            normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                        end
                        % Fill up the matrix appropriately with the sample
                        LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                            255-normFrame;
                        LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                    end
                    PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(1);
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(2);
                    FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(3);
                    oddFlipCounter
                end
            else
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesAmbigProb(gapCounter).GapID = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesOffProb(gapCounter).GapID = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCrossesOnProb(gapCounter).GapID = ...
                    0;
                oddFlipCounter
            end
            
        end
    end


    for evenFlipCounter = 1:length(FBFS.ExpNum(flyCounter).BehavData.EvenFlips)
        for gapCounter = 1:NumGaps
            if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesAbsFrameStart(gapCounter).GapID
                startFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesAbsFrameStart(gapCounter).GapID;
                endFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesAbsFrameEnd(gapCounter).GapID;
                
                for eventCounter = 1:length(startFrames)
                    startFrame  = startFrames(eventCounter);
                    endFrame    = endFrames(eventCounter);
                    
                    LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                    % Tracker for how many frames were written within sample
                    numFramesWritten = 0;
                    % Read in the appropriate frames of interest 
                    % (center frame + relativeFrames)
                    for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                        numFramesWritten = numFramesWritten + 1;
                        % Check to make sure that the frames requested are within the time points of the flip
                        % If not, fill in a white frame
                        if frameCounter < startFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        elseif frameCounter > endFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        else
                            frame = read(reader1, frameCounter);
                            % Only reads in one channel since not RGB
                            frame = frame(:,:,1);
                            % Only read in the pixels in the region of interest
                            frame = frame((GapCentsYEven(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYEven(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                         (GapCentsXEven(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXEven(gapCounter,flyCounter)+round(frameWidth/2)));
                            % Normalize the frame
                            normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                        end
                        
                        % Fill up the matrix appropriately with the sample
                        LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                            255-normFrame;
                        LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                    end
                    PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(1);
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(2);
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(3);
                    evenFlipCounter
                end
            else
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                    0;
                evenFlipCounter
            end
            
            if FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesAbsFrameStart(gapCounter).GapID
                startFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesAbsFrameStart(gapCounter).GapID;
                endFrames = FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesAbsFrameEnd(gapCounter).GapID;
                
                for eventCounter = 1:length(startFrames)
                    startFrame  = startFrames(eventCounter);
                    endFrame    = endFrames(eventCounter);
                    
                    LabeledCrossingMatInput = zeros((2*round(frameHeight/2)+1)*numFramesGrabbed,(2*round(frameWidth/2)+1));

                    % Tracker for how many frames were written within sample
                    numFramesWritten = 0;
                    % Read in the appropriate frames of interest 
                    % (center frame + relativeFrames)
                    for frameCounter = (round((startFrame+endFrame)/2) + relativeFramesGrabbed)
                        numFramesWritten = numFramesWritten + 1;
                        % Check to make sure that the frames requested are within the time points of the flip
                        % If not, fill in a white frame
                        if frameCounter < startFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        elseif frameCounter > endFrame
                            frame = 255*ones((2*round(frameHeight/2)+1),(2*round(frameWidth/2)+1));
                            normFrame = frame;
                        else
                            frame = read(reader1, frameCounter);
                            % Only reads in one channel since not RGB
                            frame = frame(:,:,1);
                            % Only read in the pixels in the region of interest
                            frame = frame((GapCentsYEven(gapCounter,flyCounter)-round(frameHeight/2)):(GapCentsYEven(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                         (GapCentsXEven(gapCounter,flyCounter)-round(frameWidth/2)):(GapCentsXEven(gapCounter,flyCounter)+round(frameWidth/2)));
                            % Normalize the frame
                            normFrame = 255*double(frame)./sum(sum(frame))*(2*round(frameHeight/2)+1)*(2*round(frameWidth/2)+1)*avgPixVal;
                        end
                        
                        % Fill up the matrix appropriately with the sample
                        LabeledCrossingMatInput(((numFramesWritten-1)*(2*round(frameHeight/2)+1)+1):((numFramesWritten)*(2*round(frameHeight/2)+1)),:) = ...
                            255-normFrame;
                        LabeledCrossingMatInput = uint8(LabeledCrossingMatInput);
                    end
                    PredAmbigOffOn = predict(netOnOffAmbig,LabeledCrossingMatInput);
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(1);
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(2);
                    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                        PredAmbigOffOn(3);
                    evenFlipCounter
                end
            else
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesAmbigProb(gapCounter).GapID(eventCounter) = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesOffProb(gapCounter).GapID(eventCounter) = ...
                    0;
                FBFS.ExpNum(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCrossesOnProb(gapCounter).GapID(eventCounter) = ...
                    0;
                evenFlipCounter
            end
            
        end
    end
end


WS.FlipBinnedFlyStruct = FBFS;

end