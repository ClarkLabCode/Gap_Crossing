% Gives the user a method to watch raw videos and use them to label a data
% set for whether or not a fly is on the glass. Can be stopped midway
% through and continued later on.

% FlipBinnedFlyStruct is the output from GapCrossingAnalysis.m
% prevStartedLabelYN should be 1 if the data is partially labeled
% LabelOrder must be accurately provided only if prevStartedLabelYN==1
% progInLabel must be accurately provided only if prevStartedLabelYN==1
% progInLabel is how many sets of flips have been done for each fly
% So if you labeled all 4 even flips and 4 odd flips for all flies in an
% experiment before pausing the script, progInLabel should be 4

function [LabeledFlyStruct, LabelOrder, progInLabel] = FlyVideoLabeler(FlipBinnedFlyStruct, prevStartedLabelYN, LabelOrder, progInLabel, LabeledFlyStruct)

% Change this to change which video you're watching. Make sure it matches
% the data in FlipBinnedFlyStruct
inputFileName = '2021-07-07 14-11-22_IsoD1_30C_berberine_chloride_40_min_training.mp4';
reader = VideoReader(inputFileName);
 
% Check if the data is already partially labeled
if prevStartedLabelYN == 1
    % If it is partially labeled, import the necessary quantities
    [NumFlies, MinFlips] = size(LabelOrder);
    NumFlies = NumFlies/2; % Dividing by two because even/odd flips
    % If labeling has started but no full set of flips is labeled,
    % reinitialize the labeled structure
    if progInLabel == 1
        LabeledFlyStruct = FlipBinnedFlyStruct;
    end
else
    % If it isn't partially labeled, initialize struct and find the quantities
    LabeledFlyStruct = FlipBinnedFlyStruct;
    progInLabel = 1;
    NumFlies = length(FlipBinnedFlyStruct.ExpNum);
    % Holds info of how many even and odd flips each fly experienced
    % Elements 1:NumFlies are number of odd flips
    % Elements NumFlies+1:2*NumFlies are number of even flips
    NumFlipsVec = zeros(1,2*NumFlies);
    for flyCounter = 1:NumFlies
        NumFlipsVec(flyCounter) = length(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips);
        NumFlipsVec(flyCounter+NumFlies) = length(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.EvenFlips);
    end
    % Find the minimum number of flips in the experiment
    MinFlips = min(NumFlipsVec);
    % Define LabelOrder which randomizes the order in which you label the
    % flips. It forces the user to label one even and odd flip for each fly
    % before moving to the next sequence.
    % Row number is which fly (first half for odd flips, second for odd)
    % Column number is the number of labeled flips so far
    LabelOrder = zeros(2*NumFlies,MinFlips);
    for flyCounter = 1:2*NumFlies
        LabelOrder(flyCounter,:) = randperm(MinFlips);
    end
end

% Now that all the quantities are found or loaded, start the main loop
for flipCounter = progInLabel:MinFlips
    for flyCounter = 1:2*NumFlies

        % Odd flips
        if flyCounter <= NumFlies
            
            % Determine the number of frames of video to watch as well as ROI
            framesInFlip = length(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).AbsoluteTime);
            minX = floor(min(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).CentroidX))-50;
            minY = floor(min(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).CentroidY))-50;
            maxX = ceil(max(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).CentroidX))+50;
            maxY = ceil(max(FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).CentroidY))+50;

            % Check that there is motion in the flip. If there isn't, then skip
            % this labeling step and just mark OffGlassYN as 0
            if minX == -50
                LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).OffGlassYN = 0;
                continue
            end

            % Initialize the parameters for playing back the video file
            startingFrameNum = FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).AbsoluteTime(1);
            defaultFrameRange = 30;     % default number of frames to play in a clip when labeling
            currentFrameRange = defaultFrameRange;
            frameNum = startingFrameNum;

            % Initialize the vector that holds whether or not a fly is off glass
            OffGlassYN = zeros(1,framesInFlip);

            % Go through every frame of the movie corresponding to the flip
            while frameNum < (framesInFlip + startingFrameNum)
                % Play back the frames from frameNum to frameNum+currentFrameRange
                for relativeFrameNum = 0:min(currentFrameRange,(framesInFlip+startingFrameNum)-frameNum)
                    frame = read(reader,frameNum + relativeFrameNum);
                    % Only show the ROI that contains the fly being labeled
                    imshow(frame(minY:maxY,minX:maxX,:));
                end
                % Request the user clicks a key to close the clip and label
                pause;
                close;
                % Ask user whether: 
                % (a) all frames in the clip had the fly on the glass
                % (b) all frames in the clip had the fly off the glass
                % (c) not all frames are the same
                % (d) the user would like the clip to be replayed
                % If an input other than (a)-(d) is selected, it defaults to replaying
                userInput = input('[a] All frames are off glass \n[b] All frames are on glass\n[c] Not all frames are the same\n[d] Replay\n','s');
                % If all frames had the fly off the glass, fill in all those
                % frames to have OffGlassYN=1
                if userInput == 'a'
                    OffGlassYN((frameNum-startingFrameNum+1):(frameNum-startingFrameNum+relativeFrameNum)) = 1;
                    % Advance the frameNum to the end of the played clip
                    frameNum = frameNum + relativeFrameNum;
                    % Default back to playing the default number of frames
                    currentFrameRange = defaultFrameRange;
                % If all frames had the fly on the glass, fill in all those
                % frames to have OffGlassYN=0
                % Note that we don't have to explicitly make OffGlassYN=0 here
                % since it was initialized to be a vector of zeros
                elseif userInput == 'b'
                    % Advance the frameNum to the end of the played clip
                    frameNum = frameNum + relativeFrameNum;
                    % Default back to playing the default number of frames
                    currentFrameRange = defaultFrameRange;
                % If not all the frames are the same, ask the user how many
                % frames they would like played at a time. This number should
                % be made repeatedly smaller until all frames in the range are
                % either all on or all off the glass
                elseif userInput == 'c'
                    currentFrameRange = input(['How many frames should we use at a time? Currently ', num2str(currentFrameRange),' frames\n']);
                end
                % By not having an outcome for any input other than (a)-(c) and
                % by not advancing the frameNum, this defaults to replaying
                
                % Check if the while loop is about to end (i.e., end of
                % flip) and give the user the option to re-label the flip
                % in case they made any errors
                if ~(frameNum < (framesInFlip + startingFrameNum))
                    relabelFlip = input('Do you want to relabel this flip? y/n\n','s');
                    % If user requests relabeling, reinitialize the
                    % parameters for this flip and restart the video playback
                    if relabelFlip == 'y'
                        % Initialize the parameters for playing back the video file
                        startingFrameNum = FlipBinnedFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).AbsoluteTime(1);
                        defaultFrameRange = 30;     % default number of frames to play in a clip when labeling
                        currentFrameRange = defaultFrameRange;
                        frameNum = startingFrameNum;

                        % Initialize the vector that holds whether or not a fly is off glass
                        OffGlassYN = zeros(1,framesInFlip);
                    end
                end
            end

            % Once the whole OffGlassYN vector is filled for a flip, add that
            % info into the LabeledFlyStruct as a new field
            LabeledFlyStruct.ExpNum(flyCounter).BehavData.OddFlips(LabelOrder(flyCounter,flipCounter)).OffGlassYN = OffGlassYN;
        
        % Even flips
        else
                        
            % Determine the number of frames of video to watch as well as ROI
            framesInFlip = length(FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).AbsoluteTime);
            minX = floor(min(FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).CentroidX))-50;
            minY = floor(min(FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).CentroidY))-50;
            maxX = ceil(max(FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).CentroidX))+50;
            maxY = ceil(max(FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).CentroidY))+50;

            % Check that there is motion in the flip. If there isn't, then skip
            % this labeling step and just mark OffGlassYN as 0
            if minX == -50
                LabeledFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).OffGlassYN = 0;
                continue
            end

            % Initialize the parameters for playing back the video file
            startingFrameNum = FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).AbsoluteTime(1);
            defaultFrameRange = 30;     % default number of frames to play in a clip when labeling
            currentFrameRange = defaultFrameRange;
            frameNum = startingFrameNum;

            % Initialize the vector that holds whether or not a fly is off glass
            OffGlassYN = zeros(1,framesInFlip);

            % Go through every frame of the movie corresponding to the flip
            while frameNum < (framesInFlip + startingFrameNum)
                % Play back the frames from frameNum to frameNum+currentFrameRange
                for relativeFrameNum = 0:min(currentFrameRange,(framesInFlip+startingFrameNum)-frameNum)
                    frame = read(reader,frameNum + relativeFrameNum);
                    % Only show the ROI that contains the fly being labeled
                    imshow(frame(minY:maxY,minX:maxX,:));
                end
                % Request the user clicks a key to close the clip and label
                pause;
                close;
                % Ask user whether: 
                % (a) all frames in the clip had the fly on the glass
                % (b) all frames in the clip had the fly off the glass
                % (c) not all frames are the same
                % (d) the user would like the clip to be replayed
                % If an input other than (a)-(d) is selected, it defaults to replaying
                userInput = input('[a] All frames are off glass \n[b] All frames are on glass\n[c] Not all frames are the same\n[d] Replay\n','s');
                % If all frames had the fly off the glass, fill in all those
                % frames to have OffGlassYN=1
                if userInput == 'a'
                    OffGlassYN((frameNum-startingFrameNum+1):(frameNum-startingFrameNum+relativeFrameNum)) = 1;
                    % Advance the frameNum to the end of the played clip
                    frameNum = frameNum + relativeFrameNum;
                    % Default back to playing the default number of frames
                    currentFrameRange = defaultFrameRange;
                % If all frames had the fly on the glass, fill in all those
                % frames to have OffGlassYN=0
                % Note that we don't have to explicitly make OffGlassYN=0 here
                % since it was initialized to be a vector of zeros
                elseif userInput == 'b'
                    % Advance the frameNum to the end of the played clip
                    frameNum = frameNum + relativeFrameNum;
                    % Default back to playing the default number of frames
                    currentFrameRange = defaultFrameRange;
                % If not all the frames are the same, ask the user how many
                % frames they would like played at a time. This number should
                % be made repeatedly smaller until all frames in the range are
                % either all on or all off the glass
                elseif userInput == 'c'
                    currentFrameRange = input(['How many frames should we use at a time? Currently ', num2str(currentFrameRange),' frames\n']);
                end
                % By not having an outcome for any input other than (a)-(c) and
                % by not advancing the frameNum, this defaults to replaying
                
                % Check if the while loop is about to end (i.e., end of
                % flip) and give the user the option to re-label the flip
                % in case they made any errors
                if ~(frameNum < (framesInFlip + startingFrameNum))
                    relabelFlip = input('Do you want to relabel this flip? y/n\n','s');
                    % If user requests relabeling, reinitialize the
                    % parameters for this flip and restart the video playback
                    if relabelFlip == 'y'
                        % Initialize the parameters for playing back the video file
                        startingFrameNum = FlipBinnedFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).AbsoluteTime(1);
                        defaultFrameRange = 30;     % default number of frames to play in a clip when labeling
                        currentFrameRange = defaultFrameRange;
                        frameNum = startingFrameNum;

                        % Initialize the vector that holds whether or not a fly is off glass
                        OffGlassYN = zeros(1,framesInFlip);
                    end
                end
            end

            % Once the whole OffGlassYN vector is filled for a flip, add that
            % info into the LabeledFlyStruct as a new field
            LabeledFlyStruct.ExpNum(flyCounter-NumFlies).BehavData.EvenFlips(LabelOrder(flyCounter,flipCounter)).OffGlassYN = OffGlassYN;
        end
    end
    
    % Increase the counter on progInLabel since one full set of flips has
    % been labeled now
    progInLabel = progInLabel+1;
    
    % Give the user the option to stop labeling for this session
    pauseScript = input('Do you want to take a break from labeling and save your progress? y/n\n','s');
    if pauseScript == 'y'
        break
    end
    
end

end