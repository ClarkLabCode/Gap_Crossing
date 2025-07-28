% VISUALIZERAWDATASCRIPTS Script to generate many visualizations of raw data
%  
%  Currently plots things such as trajectory in the cassette, distribution
%  of velocities as a function of y position, angular velocities vs y, etc.
%
%  Some of the scripts do this for an individual fly, while others do it by
%  randomly subsampling data from all flies.
% 
%  This script can be used as a template for plotting other visualizations
%  of the raw data.

%% Setting up stuff for plotting examples later

% Grab FBFS from WS
FBFS = WS.FlipBinnedFlyStruct;

% All the skeletons look basically the same after alignment, so just take
% the first one
SkelX = FBFS(1).IdData.AlignedSkeleton_x;
SkelY = FBFS(1).IdData.AlignedSkeleton_y;

% Initialize some variables
maxFlips = 0; % Variable for tracking total flips in experiment
numFlies = length(FBFS);
NumGaps = 4; % Everything written here only currently works for 4 gaps

% Go through and see how many flips were in the experiments by taking the
% largest number of flips from all flies
for flyCounter = 1:numFlies
    maxFlips = max([maxFlips, length(FBFS(flyCounter).AlignedData.OddFlips), length(FBFS(flyCounter).AlignedData.EvenFlips)]);
end

% Set rng seed so we can replicate things later
s = rng(0);
flipsPerFly = 10; % Number of flips we'll plot from each fly
randFlips = randi(maxFlips,[flipsPerFly,numFlies]); % Matrix of random flips we'll plot from each fly
minFramesInFlip = 228; % Minimum number of frames a fly must be tracked for within a flip to be plotted

% Initialize some matrices that will hold all of the data we have
% 4D matrices with following structure: (frames in flip, flip #, fly #, odd/even flip)
IsoD1MatrixTime = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixX = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixY = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixTheta = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixMajorAxisLength = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixMinorAxisLength = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixVelX = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixVelY = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixVelTheta = NaN(minFramesInFlip,maxFlips,numFlies,2);

% Now fill in the matrices we initialized above with all the data we have
for flyCounter = 1:numFlies
    % Do it for odd flips
    for oddFlipCounter = 1:length(FBFS(flyCounter).AlignedData.OddFlips)
        % We're only plotting the data from frames with at least as many
        % frames as minFramesInFlip, so don't change the NaN if it doesn't
        % meet that criteria
        if length(FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).FrameInFlip) >= minFramesInFlip+1
            % Grab all the locomotion data
            IsoD1MatrixTime(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip(1:minFramesInFlip);
            IsoD1MatrixX(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedCentroidX(1:minFramesInFlip);
            IsoD1MatrixY(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidY(1:minFramesInFlip);
            IsoD1MatrixTheta(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVerticalTheta(1:minFramesInFlip);
            IsoD1MatrixMajorAxisLength(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).MajorAxisLength(1:minFramesInFlip);
            IsoD1MatrixMinorAxisLength(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).MinorAxisLength(1:minFramesInFlip);
            IsoD1MatrixVelX(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelX(1:minFramesInFlip);
            IsoD1MatrixVelY(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelY(1:minFramesInFlip);
            IsoD1MatrixVelTheta(:,oddFlipCounter,flyCounter,1) = ...
                FBFS(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelTheta(1:minFramesInFlip);
        end
    end
    % Do it for even flips
    for evenFlipCounter = 1:length(FBFS(flyCounter).AlignedData.EvenFlips)
        % We're only plotting the data from frames with at least as many
        % frames as minFramesInFlip, so don't change the NaN if it doesn't
        % meet that criteria
        if length(FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).FrameInFlip) >= minFramesInFlip+1
            IsoD1MatrixTime(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedTimeInFlip(1:minFramesInFlip);
            IsoD1MatrixX(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedCentroidX(1:minFramesInFlip);
            IsoD1MatrixY(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidY(1:minFramesInFlip);
            IsoD1MatrixTheta(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVerticalTheta(1:minFramesInFlip);
            IsoD1MatrixMajorAxisLength(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).MajorAxisLength(1:minFramesInFlip);
            IsoD1MatrixMinorAxisLength(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).MinorAxisLength(1:minFramesInFlip);
            IsoD1MatrixVelX(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelX(1:minFramesInFlip);
            IsoD1MatrixVelY(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelY(1:minFramesInFlip);
            IsoD1MatrixVelTheta(:,evenFlipCounter,flyCounter,2) = ...
                FBFS(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelTheta(1:minFramesInFlip);
        end
    end
end

% Now initialize a bunch of things to use for identifying what behavior
% took place at each of the gap events. Much of this code is copied from
% the algorithm and method used in FindCrossEvents.m
% This part may take many read throughs to fully follow/remember what's
% being done, sorry about that... I've done my best to write out all the
% details in the comments but am happy to answer questions about this at a
% later time if the comments are insufficient (from Joe)

% Initialize 6D behavior matrix with following structure:
% (Start/End Frame #, flip #, fly #, odd/even flip, gap #, behavior #)
% Behavior #: 1 = Proper Crossing, 2 = Glass Crossing, 3 = Circum, 4 = Retreat
% Note that the first dimension is size 10*2 to allow for up to 10 occurrences 
% of any singular behavior within any flip. The structuring within the 
% 20 is [start occ1, end occ1, start occ2, end occ2, ..., start occ10, end occ10]
IsoD1MatrixBehav = NaN(10*2,maxFlips,numFlies,2,NumGaps,4);

% Make the vectors that hold the appropriate compartment numbers for each
% type of region of interest (see corridor diagram for compartment labels)
Gaps = 2*(1:NumGaps);
Corrs = 2*(0:NumGaps)+1;
Wells = (2*NumGaps + 1) + (1:NumGaps);

% 1st dimension is CompID, 2nd is GapNumber, 3rd is For/Back
% Forward is defined as along the direction from smallest to largest gap
CrossID = zeros(3,NumGaps,2);
RetID = zeros(3,NumGaps,2);
CircID = zeros(3,NumGaps,2);

% Define the transition IDs for each event of interest (1 = For, 2 = Back)
% Forward is defined as along the direction from smallest to largest gap
for GapNumber = 1:NumGaps
    CrossID(:,GapNumber,1) = [Corrs(GapNumber), Gaps(GapNumber), Corrs(GapNumber+1)];
    CrossID(:,GapNumber,2) = [Corrs(GapNumber+1), Gaps(GapNumber), Corrs(GapNumber)];
    RetID(:,GapNumber,1) = [Corrs(GapNumber), Gaps(GapNumber), Corrs(GapNumber)];
    RetID(:,GapNumber,2) = [Corrs(GapNumber+1), Gaps(GapNumber), Corrs(GapNumber+1)];
    CircID(:,GapNumber,1) = [Corrs(GapNumber), Gaps(GapNumber), Wells(GapNumber)];
    CircID(:,GapNumber,2) = [Corrs(GapNumber+1), Gaps(GapNumber), Wells(GapNumber)];
end

trackerVar = 0;
bufferFrames = 0;

% Go through every flip of every fly and grab the start and end frames of
% each gap event by looking through UniqCompID and UniqCompIDIndex with the
% help of the above defined transition IDs
for oddEven = 1:2
    for flyCounter = 1:numFlies
        if oddEven == 1
            % When doing the odd loop, grab the OddFlips data
            FlipsData = FBFS(flyCounter).AlignedData.OddFlips;
        else
            % When doing the even loop, grab the EvenFlips data
            FlipsData = FBFS(flyCounter).AlignedData.EvenFlips;
        end
        for flipCounter = 1:length(FlipsData)
            for GapNumber = 1:NumGaps
                % Continually update UniqCompIDIndex and UniqCompID
                % throughout the loop so that we don't have to dot index as
                % many things on each line
                UniqCompIDIndex = FlipsData(flipCounter).UniqCompIDIndex;
                UniqCompID = FlipsData(flipCounter).UniqCompID;
                % Grab the frame(s) during which cross event(s) starts
                % This is found by computing the last frame in which the
                % fly leaves the corridor and enters the gap (hence the +1
                % within the index of UniqCompIDIndex followed by the -1
                % outside of it). 
                % Note that this may return several frames because there
                % may be several occurrences of a cross event within a
                % singular flip.
                startFrame = ...
                    UniqCompIDIndex(strfind(UniqCompID, CrossID(:,GapNumber,oddEven)')+1)-1;
                % Trying to shift the start frame to be bufferFrames frames earlier
                if length(startFrame) > 0
                    startFrame = startFrame-bufferFrames;
                    startFrame(startFrame<1) = 1;
                end
                % Grab the frame(s) during which cross event(s) ends
                % This is found by computing the first frame in which the
                % fly leaves the gap and enters the corridor (hence the +2
                % within the index of UniqCompIDIndex).
                % Again, note that this may return several frames because
                % there may be several occurrences of a cross event within a
                % singular flip.
                endFrame = ...
                    UniqCompIDIndex(strfind(UniqCompID, CrossID(:,GapNumber,oddEven)')+2);
                % Assuming there are any occurences at all of crossing
                % events within the flip, go through each occurrence
                if length(startFrame) > 0
                    % endFrame must always be at most the last frame we are
                    % grabbing from each flip, so use the min function to
                    % ensure this.
                    % Note that this must be done within the if statement
                    % because min will override the NaNs otherwise.
                    endFrame = min(endFrame+bufferFrames,minFramesInFlip);
                    for occurrenceCounter = 1:length(startFrame)
                        % There's a further complication when grabbing cross events
                        % because we must determine whether they were a proper or a
                        % glass cross. To do this, check the renormalized NN output
                        % and assign the event as a proper or glass cross accordingly
                        % If NN output is < 0.5, it's a proper crossing
                        if FlipsData(flipCounter).RenormUpGlassCrossProb(GapNumber).GapID(occurrenceCounter) < 0.5
                            % Fill in start frame(s) for PC
                            IsoD1MatrixBehav(2*(occurrenceCounter-1)+1,flipCounter,flyCounter,oddEven,GapNumber,1) = ...
                                startFrame(occurrenceCounter);
                            % Fill in end frame(s) for PC
                            IsoD1MatrixBehav(2*(occurrenceCounter-1)+2,flipCounter,flyCounter,oddEven,GapNumber,1) = ...
                                endFrame(occurrenceCounter);
                        % If NN output is >= 0.5, it's a glass crossing
                        else
                            % Fill in start frame(s) for GC
                            IsoD1MatrixBehav(2*(occurrenceCounter-1)+1,flipCounter,flyCounter,oddEven,GapNumber,2) = ...
                                startFrame(occurrenceCounter);
                            % Fill in end frame(s) for GC
                            IsoD1MatrixBehav(2*(occurrenceCounter-1)+2,flipCounter,flyCounter,oddEven,GapNumber,2) = ...
                                endFrame(occurrenceCounter);
                        end
                    end
                end
                % Grab the frame(s) during which circum event(s) starts
                % This is found by computing the last frame in which the
                % fly leaves the corridor and enters the gap (hence the +1
                % within the index of UniqCompIDIndex followed by the -1
                % outside of it). 
                % Note that this may return several frames because there
                % may be several occurrences of a circum event within a
                % singular flip.
                startFrame = ...
                    UniqCompIDIndex(strfind(UniqCompID, CircID(:,GapNumber,oddEven)')+1)-1;
                % Trying to shift the start frame to be bufferFrames frames earlier
                if length(startFrame) > 0
                    startFrame = startFrame-bufferFrames;
                    startFrame(startFrame<1) = 1;
                end
                % Grab the frame(s) during which circum event(s) ends
                % This is found by computing the first frame in which the
                % fly leaves the gap for the second time (hence the +4
                % within the index of UniqCompIDIndex).
                % Again, note that this may return several frames because
                % there may be several occurrences of a circum event within a
                % singular flip.
                % As one further subtlty, note the need for the min
                % function within UniqCompIDIndex because it is possible
                % for the final frame of a circumvention to take place in a
                % frame beyond the final frame of UniqCompID, so in those
                % (rare) cases, we just take the final frame from
                % UniqCompID to be our final frame rather than generate an
                % error.
                endFrame = ...
                    UniqCompIDIndex(min(strfind(UniqCompID, CircID(:,GapNumber,oddEven)')+4,length(UniqCompID)));
                % Assuming there are any occurences at all of circum
                % events within the flip, go through each occurrence
                if length(startFrame) > 0
                    % endFrame must always be at most the last frame we are
                    % grabbing from each flip, so use the min function to
                    % ensure this.
                    % Note that this must be done within the if statement
                    % because min will override the NaNs otherwise.
                    endFrame = min(endFrame+bufferFrames,minFramesInFlip);
                    for occurrenceCounter = 1:length(startFrame)
                        % Fill in start frame(s) for Circ
                        IsoD1MatrixBehav(2*(occurrenceCounter-1)+1,flipCounter,flyCounter,oddEven,GapNumber,3) = ...
                            startFrame(occurrenceCounter);
                        % Fill in end frame(s) for Circ
                        IsoD1MatrixBehav(2*(occurrenceCounter-1)+2,flipCounter,flyCounter,oddEven,GapNumber,3) = ...
                            endFrame(occurrenceCounter);
                    end
                end
                % Grab the frame(s) during which ret event(s) starts
                % This is found by computing the last frame in which the
                % fly leaves the corridor and enters the gap (hence the +1
                % within the index of UniqCompIDIndex followed by the -1
                % outside of it). 
                % Note that this may return several frames because there
                % may be several occurrences of a ret event within a
                % singular flip.
                startFrame = ...
                    UniqCompIDIndex(strfind(UniqCompID, RetID(:,GapNumber,oddEven)')+1)-1;
                % Trying to shift the start frame to be bufferFrames frames earlier
                if length(startFrame) > 0
                    startFrame = startFrame-bufferFrames;
                    startFrame(startFrame<1) = 1;
                end
                % Grab the frame(s) during which ret event(s) ends
                % This is found by computing the first frame in which the
                % fly leaves the gap and enters the corridor (hence the +2
                % within the index of UniqCompIDIndex).
                % Again, note that this may return several frames because
                % there may be several occurrences of a ret event within a
                % singular flip.
                endFrame = ...
                    UniqCompIDIndex(strfind(UniqCompID, RetID(:,GapNumber,oddEven)')+2);
                % Assuming there are any occurences at all of ret
                % events within the flip, go through each occurrence
                if length(startFrame) > 0
                    % endFrame must always be at most the last frame we are
                    % grabbing from each flip, so use the min function to
                    % ensure this.
                    % Note that this must be done within the if statement
                    % because min will override the NaNs otherwise.
                    endFrame = min(endFrame+bufferFrames,minFramesInFlip);
                    for occurrenceCounter = 1:length(startFrame)
                        % Fill in start frame(s) for Ret
                        IsoD1MatrixBehav(2*(occurrenceCounter-1)+1,flipCounter,flyCounter,oddEven,GapNumber,4) = ...
                            startFrame(occurrenceCounter);
                        % Fill in end frame(s) for Ret
                        IsoD1MatrixBehav(2*(occurrenceCounter-1)+2,flipCounter,flyCounter,oddEven,GapNumber,4) = ...
                            endFrame(occurrenceCounter);
                    end
                end
            end
        end
    end
end

% Construct a straightened skeleton from the raw skeletons
StraightenedSkelX = zeros(size(SkelX));
StraightenedSkelY = zeros(size(SkelY));

StraightenedSkelX([1,2,5,6,9,10,13,14,17,18],:) = repmat(min(SkelX([1,2,5,6,9,10,13,14,17,18],:)),size([1,2,5,6,9,10,13,14,17,18]'));
StraightenedSkelX([3,4,7,8,11,12,15,16],:) = repmat(min(SkelX([3,4,7,8,11,12,15,16],:)),size([3,4,7,8,11,12,15,16]'));
StraightenedSkelX(18+[1,2,5,6,9,10,13,14,17,18],:) = repmat(max(SkelX(18+[1,2,5,6,9,10,13,14,17,18],:))+4/15,size([1,2,5,6,9,10,13,14,17,18]'));
StraightenedSkelX(18+[3,4,7,8,11,12,15,16],:) = repmat(max(SkelX(18+[2,5,6,9,10,13,14,17],:))+5+4/15,size([3,4,7,8,11,12,15,16]'));

StraightenedSkelY([1,36],:) = repmat(min(SkelY([1,36],:)),size([1,36]'));
StraightenedSkelY([18,19],:) = repmat(max(SkelY([18,19],:)),size([18,19]'));
StraightenedSkelY([2,3,34,35],:) = repmat(min(SkelY([2,3,34,35],:)),size([2,3,34,35]'));
StraightenedSkelY([4,5,32,33],:) = repmat(max(SkelY([4,5,32,33],:)),size([4,5,32,33]'));
StraightenedSkelY([6,7,30,31],:) = repmat(min(SkelY([6,7,30,31],:)),size([6,7,30,31]'));
StraightenedSkelY([8,9,28,29],:) = repmat(max(SkelY([8,9,28,29],:)),size([8,9,28,29]'));
StraightenedSkelY([10,11,26,27],:) = repmat(min(SkelY([10,11,26,27],:)),size([10,11,26,27]'));
StraightenedSkelY([12,13,24,25],:) = repmat(max(SkelY([12,13,24,25],:)),size([12,13,24,25]'));
StraightenedSkelY([14,15,22,23],:) = repmat(min(SkelY([14,15,22,23],:)),size([14,15,22,23]'));
StraightenedSkelY([16,17,20,21],:) = repmat(max(SkelY([16,17,20,21],:)),size([16,17,20,21]'));


%% Plot all Y vs T of a given fly

% Example of how to plot all odd/even flips of a given fly
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
chosenFlyNum = 1;
oddEven = 1;
figure
plot(IsoD1MatrixTime(:,:,chosenFlyNum,oddEven),...
     IsoD1MatrixY(:,:,chosenFlyNum,oddEven),...
     'Color',[0.7,0.7,0.7])
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
hold off

title('y(t) for all flips of one fly')
xlabel('Time in flip (sec)')
ylabel('Y position (mm)')

%% Plot random subset of Y vs T from all flies

% Example of how to plot the randomly selected odd/even flips across all flies
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
oddEven = 1;
figure
hold on
for flyCounter = 1:numFlies
    plot(IsoD1MatrixTime(:,randFlips(:,flyCounter),flyCounter,oddEven),...
         IsoD1MatrixY(:,randFlips(:,flyCounter),flyCounter,oddEven),...
         'Color',[0.7,0.7,0.7])
end
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
hold off

title('y(t) for random flips of all flies')
xlabel('Time in flip (sec)')
ylabel('Y position (mm)')

%% Plot Y vs T for flips starting below a certain point on corridor

% If we want to only plot flips in which a fly starts below a certain
% position, set this threshold position and find flips in which this happens
lowPositionThresh = -16.5; % The position below which we plot a flip
% Initialize a logical array that'll be used to filter what to plot
startsLow = zeros(maxFlips,numFlies,2);
fracOfLowStarts = 1; % Sets the fraction of flips to display that meet the criteria
% Go through each flip of each fly and determine whether or not it starts
% below the position threshold. If so, fill startsLow with a 1.
for flyCounter = 1:numFlies
    for flipCounter = 1:maxFlips
        for oddEvenCounter = 1:2
            if IsoD1MatrixY(1,flipCounter,flyCounter,oddEvenCounter) < lowPositionThresh
                % To implement fracOfLowStarts, use a rng to determine whether
                % or not it'll be plotted
                if rand(1) < fracOfLowStarts
                    startsLow(flipCounter,flyCounter,oddEvenCounter) = 1;
                end
            end
        end
    end
end
% Convert startsLow to logical matrices
startsLow = logical(startsLow);

% Example of how to plot odd/even flips that start below a given position
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
oddEven = 1;
figure
% This is the version of the logical array that will have the flip parity
% that we aren't interested in zeroed out so that we don't plot those
startsLowOddEven = startsLow;
% The 3-oddEven switches the odd/even logic (i.e., when we want to plot odd
% flips which corresponds to 1, we should zero out all the even flips and
% vice versa)
startsLowOddEven(:,:,3-oddEven) = 0;
hold on
plot(IsoD1MatrixTime(:,startsLowOddEven),...
         IsoD1MatrixY(:,startsLowOddEven),...
         'Color',[0.7,0.7,0.7])
% Overlay the median trajectory in black
plot(nanmedian(IsoD1MatrixTime(:,startsLowOddEven),2),...
     nanmedian(IsoD1MatrixY(:,startsLowOddEven),2),'k')

% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
hold off

title('y(t) for some flips starting below certain y position of all flies')
xlabel('Time in flip (sec)')
ylabel('Y position (mm)')

%% Plot Y vs Theta for a given fly
% Example of how to plot all odd/even flips of a given fly
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
chosenFlyNum = 1;
oddEven = 1;
figure
plot(IsoD1MatrixTheta(:,:,chosenFlyNum,oddEven),...
     IsoD1MatrixY(:,:,chosenFlyNum,oddEven),...
     'Color',[0.7,0.7,0.7])
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--')
hold off

title('$\theta(y)$ for all flips of one fly','Interpreter','latex')
xlabel('Orientation (degrees)')
ylabel('Y position (mm)')

%% Plot Y vs VelTheta for a given fly
% Example of how to plot all odd/even flips of a given fly
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
chosenFlyNum = 1;
oddEven = 1;
figure
plot(IsoD1MatrixVelTheta(:,:,chosenFlyNum,oddEven),...
     IsoD1MatrixY(:,:,chosenFlyNum,oddEven),...
     'Color',[0.7,0.7,0.7])
% Example of how to overlay the gap locations in the plots
xl = xlim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--')
hold off

title('$\omega(y)$ for all flips of one fly','Interpreter','latex')
xlabel('Angular speed (degrees/sec)')
ylabel('Y position (mm)')

%% Plot Y vs VelY for a given fly
% Example of how to plot all odd/even flips of a given fly
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
chosenFlyNum = 1;
oddEven = 1;
figure
plot(IsoD1MatrixVelY(:,:,chosenFlyNum,oddEven),...
     IsoD1MatrixY(:,:,chosenFlyNum,oddEven),...
     'Color',[0.7,0.7,0.7])
% Example of how to overlay the gap locations in the plots
xl = xlim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--')
hold off

title('$\dot{y}(y)$ for all flips of one fly','Interpreter','latex')
xlabel('Vertical velocity (mm/sec)')
ylabel('Y position (mm)')
xlim([-50,50])

%% Plot Y vs X for a given fly
% Example of how to plot all odd/even flips of a given fly
% For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
chosenFlyNum = 3;
oddEven = 1;
figure
plot(IsoD1MatrixX(:,:,chosenFlyNum,oddEven),...
     IsoD1MatrixY(:,:,chosenFlyNum,oddEven),...
     'Color',[0.7,0.7,0.7])
pbaspect([1080 1920 1])
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--')
hold off

title('Trajectory of all flips of one fly','Interpreter','latex')
xlabel('X position (mm)')
ylabel('Y position (mm)')

% Plot Y vs VelY in a histogram for all flies
figure
numBins = 200;
oddEven = 1;
histX = reshape([IsoD1MatrixVelY(:,:,:,oddEven)],[],1);
histX = histX(~isnan(histX));
histY = reshape([IsoD1MatrixY(:,:,:,oddEven)],[],1);
histY = histY(~isnan(histY));
XBinLim = 40;
YBinLim = 24;
histogram2(histX,histY,-XBinLim:(2*XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
set(gca,'ColorScale','log');
colorbar;
grid off
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    a = plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--');
    uistack(a,'top');
end
a = plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--');
uistack(a,'top');
hold off
xlabel('Vertical Velocity (mm/sec)');
ylabel('Y Position (mm)');

% Plot Y vs Theta in a histogram for all flies
figure
numBins = 200;
oddEven = 1;
histX = reshape([IsoD1MatrixTheta(:,:,:,oddEven)],[],1);
histX = histX(~isnan(histX));
histY = reshape([IsoD1MatrixY(:,:,:,oddEven)],[],1);
histY = histY(~isnan(histY));
XBinLim = 100;
YBinLim = 24;
histogram2(histX,histY,-XBinLim:(2*XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
set(gca,'ColorScale','log');
colorbar;
grid off
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    a = plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--');
    uistack(a,'top');
end
a = plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--');
uistack(a,'top');
hold off
xlabel('Orientation (degrees)');
ylabel('Y Position (mm)');

% Plot Y vs VelTheta in a histogram for all flies
figure
numBins = 200;
oddEven = 1;
histX = reshape([IsoD1MatrixVelTheta(:,:,:,oddEven)],[],1);
histX = histX(~isnan(histX));
histY = reshape([IsoD1MatrixY(:,:,:,oddEven)],[],1);
histY = histY(~isnan(histY));
XBinLim = 2500;
YBinLim = 24;
histogram2(histX,histY,-XBinLim:(2*XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
set(gca,'ColorScale','log');
colorbar;
grid off
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    a = plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--');
    uistack(a,'top');
end
a = plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--');
uistack(a,'top');
hold off
xlabel('Angular Speed (degrees/sec)');
ylabel('Y Position (mm)');

%% Doing a bunch of stuff down here with filtering by specific behaviors
dimSizes = size(IsoD1MatrixX);
ReshapedIsoD1MatrixX = zeros(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedIsoD1MatrixY = zeros(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedIsoD1MatrixVelY = zeros(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedIsoD1MatrixTime = zeros(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedIsoD1MatrixVelTheta = zeros(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
for i = 1:dimSizes(1)
    for k = 1:dimSizes(3)
        for j = 1:dimSizes(2)
            for m = 1:dimSizes(4)
                ReshapedIsoD1MatrixX(i,(k-1)*dimSizes(2)+j,m) = ...
                    IsoD1MatrixX(i,j,k,m);
                ReshapedIsoD1MatrixY(i,(k-1)*dimSizes(2)+j,m) = ...
                    IsoD1MatrixY(i,j,k,m);
                ReshapedIsoD1MatrixVelY(i,(k-1)*dimSizes(2)+j,m) = ...
                    IsoD1MatrixVelY(i,j,k,m);
                ReshapedIsoD1MatrixTime(i,(k-1)*dimSizes(2)+j,m) = ...
                    IsoD1MatrixTime(i,j,k,m);
                ReshapedIsoD1MatrixVelTheta(i,(k-1)*dimSizes(2)+j,m) = ...
                    IsoD1MatrixVelTheta(i,j,k,m);
            end
        end
    end
end

dimSizesBehav = size(IsoD1MatrixBehav);
ReshapedIsoD1MatrixBehav = ...
    zeros(dimSizesBehav(1),dimSizesBehav(2)*dimSizesBehav(3),dimSizesBehav(4),dimSizesBehav(5),dimSizesBehav(6));
for i = 1:dimSizesBehav(1)
    for k = 1:dimSizesBehav(3)
        for j = 1:dimSizesBehav(2)
            for m = 1:dimSizesBehav(4)
                for n = 1:dimSizesBehav(5)
                    for o = 1:dimSizesBehav(6)
                        ReshapedIsoD1MatrixBehav(i,(k-1)*dimSizesBehav(2)+j,m,n,o) = ...
                            IsoD1MatrixBehav(i,j,k,m,n,o);
                    end
                end
            end
        end
    end
end


for BehavNum = 1:4

% figure


IsoD1MatrixIndexing = false(size(IsoD1MatrixX));

for chosenFlyNum = 1:28

for GapNum = 1:4


    for flyCounter = 1:numFlies
        for flipCounter = 1:maxFlips
            if ~isnan(IsoD1MatrixBehav(1,flipCounter,flyCounter,oddEven,GapNum,BehavNum))
                IsoD1MatrixIndexing(IsoD1MatrixBehav(1,flipCounter,flyCounter,oddEven,GapNum,BehavNum):...
                                    IsoD1MatrixBehav(2,flipCounter,flyCounter,oddEven,GapNum,BehavNum),...
                                    flipCounter,flyCounter,oddEven) = 1;
            end
        end
    end

end
end

    FilteredIsoD1Time = IsoD1MatrixTime.*IsoD1MatrixIndexing;
    FilteredIsoD1TimePerm = FilteredIsoD1Time;
    FilteredIsoD1Time(FilteredIsoD1TimePerm == 0) = NaN;
    FilteredIsoD1X = IsoD1MatrixX.*IsoD1MatrixIndexing;
    FilteredIsoD1X(FilteredIsoD1TimePerm == 0) = NaN;
    FilteredIsoD1Y = IsoD1MatrixY.*IsoD1MatrixIndexing;
    FilteredIsoD1Y(FilteredIsoD1TimePerm == 0) = NaN;
    FilteredIsoD1Theta = IsoD1MatrixTheta.*IsoD1MatrixIndexing;
    FilteredIsoD1Theta(FilteredIsoD1TimePerm == 0) = NaN;
    FilteredIsoD1VelX = IsoD1MatrixVelX.*IsoD1MatrixIndexing;
    FilteredIsoD1VelX(FilteredIsoD1TimePerm == 0) = NaN;
    FilteredIsoD1VelY = IsoD1MatrixVelY.*IsoD1MatrixIndexing;
    FilteredIsoD1VelY(FilteredIsoD1TimePerm == 0) = NaN;
    FilteredIsoD1VelTheta = IsoD1MatrixVelTheta.*IsoD1MatrixIndexing;
    FilteredIsoD1VelTheta(FilteredIsoD1TimePerm == 0) = NaN;


ReshapedFilteredIsoD1Time = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedFilteredIsoD1X = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedFilteredIsoD1Y = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedFilteredIsoD1Theta = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedFilteredIsoD1VelX = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedFilteredIsoD1VelY = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));
ReshapedFilteredIsoD1VelTheta = NaN(dimSizes(1),dimSizes(2)*dimSizes(3),dimSizes(4));

for i = 1:dimSizes(1)
    for k = 1:dimSizes(3)
        for j = 1:dimSizes(2)
            for m = 1:dimSizes(4)
                ReshapedFilteredIsoD1Time(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1Time(i,j,k,m);
                ReshapedFilteredIsoD1X(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1X(i,j,k,m);
                ReshapedFilteredIsoD1Y(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1Y(i,j,k,m);
                ReshapedFilteredIsoD1Theta(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1Theta(i,j,k,m);
                ReshapedFilteredIsoD1VelX(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1VelX(i,j,k,m);
                ReshapedFilteredIsoD1VelY(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1VelY(i,j,k,m);
                ReshapedFilteredIsoD1VelTheta(i,(k-1)*dimSizes(2)+j,m) = ...
                    FilteredIsoD1VelTheta(i,j,k,m);
            end
        end
    end
end

    % IsoD1MatrixIndexing = false(size(ReshapedIsoD1MatrixX));
    % 
    % for counter = 1:(dimSizesBehav(2)*dimSizesBehav(3))
    %     if ~isnan(ReshapedIsoD1MatrixBehav(1,counter,oddEven,GapNum,BehavNum))
    %         IsoD1MatrixIndexing(ReshapedIsoD1MatrixBehav(1,counter,oddEven,GapNum,BehavNum):...
    %                             ReshapedIsoD1MatrixBehav(2,counter,oddEven,GapNum,BehavNum),...
    %                             counter,oddEven) = 1;
    %     end
    % end
    % 
    % FilteredIsoD1X = ReshapedIsoD1MatrixX.*IsoD1MatrixIndexing;
    % FilteredIsoD1X(FilteredIsoD1X == 0) = NaN;
    % FilteredIsoD1Y = ReshapedIsoD1MatrixY.*IsoD1MatrixIndexing;
    % FilteredIsoD1Y(FilteredIsoD1Y == 0) = NaN;

% for chosenFlyNum = 1:28
%     hold on
%     % Example of how to plot all odd/even flips of a given fly
%     % For odd flips, set oddEven to 1. For even flips, set oddEven to 2.
%     % chosenFlyNum = 1;
%     oddEven = 1;
%     plot(FilteredIsoD1X(:,:,chosenFlyNum,oddEven),...
%     FilteredIsoD1Y(:,:,chosenFlyNum,oddEven),...
%     'Color',[0.7,0.7,0.7])
%     pbaspect([1080 1920 1])
%     % Example of how to overlay the gap locations in the plots
%     xl = [0,8];
%     yl = ylim;
%     for i = [1,2,5,6,9,10,13,14,17,18]
%         plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
%     end
%     plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--')
%     title('Trajectory of all flips of one fly','Interpreter','latex')
%     xlabel('X position (mm)')
%     ylabel('Y position (mm)')
%     hold off
% end

% Plot Y vs VelY in a histogram for all flies
figure
numBins = 200;
oddEven = 1;
histX = reshape([FilteredIsoD1VelTheta(:,:,:,oddEven)],[],1);
histX = histX(~isnan(histX));
histY = reshape([FilteredIsoD1Y(:,:,:,oddEven)],[],1);
histY = histY(~isnan(histY));
XBinLim = 2500;
YBinLim = 24;
histogram2(histX,histY,-XBinLim:(2*XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
% histogram2(histX,histY,0:(XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
set(gca,'ColorScale','log');
% pbaspect([1080 1920 1])
colorbar;
grid off
% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
for i = [1,2,5,6,9,10,13,14,17,18]
    a = plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--');
    uistack(a,'top');
end
a = plot(zeros(11,1),yl(1):(yl(2)-yl(1))/10:yl(2),'k--');
uistack(a,'top');
hold off
xlabel('Angular Velocity (deg/sec)');
ylabel('Y Position (mm)');
switch BehavNum
    case 1
        title('Proper Crossing')
    case 2
        title('Glass Crossing')
    case 3
        title('Circumvention')
    case 4
        title('Retreat')
end

% Plotting trajectory for all different gap events across all flies

figure

% numSamples = size(ReshapedFilteredIsoD1X,2);
numSamples = 100;
randSubset = randsample(size(ReshapedFilteredIsoD1X,2),numSamples);
plot(ReshapedFilteredIsoD1X(:,randSubset,1),ReshapedFilteredIsoD1Y(:,randSubset,1))
xlim([0,8])
% Example of how to overlay the skeleton in the plots
xl = xlim;
yl = ylim;
hold on
a = plot([0;StraightenedSkelX(end/2+1:end,1);0],...
         [StraightenedSkelY(19,1);StraightenedSkelY(end/2+1:end,1);StraightenedSkelY(36,1)],'k');
uistack(a,'top');
hold off
pbaspect([1080 1920 1])
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
switch BehavNum
    case 1
        title('Proper Crossing')
    case 2
        title('Glass Crossing')
    case 3
        title('Circumvention')
    case 4
        title('Retreat')
end

end