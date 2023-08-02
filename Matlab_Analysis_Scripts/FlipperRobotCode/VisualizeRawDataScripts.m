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

% All the skeletons look basically the same after alignment, so just take
% the first one
SkelX = WS.FlipBinnedFlyStruct(1).IdData.AlignedSkeleton_x;
SkelY = WS.FlipBinnedFlyStruct(1).IdData.AlignedSkeleton_y;

% Initialize some variables
maxFlips = 0; % Variable for tracking total flips in experiment
numFlies = length(WS.FlipBinnedFlyStruct);

% Go through and see how many flips were in the experiments by taking the
% largest number of flips from all flies
for flyCounter = 1:numFlies
    maxFlips = max([maxFlips, length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips), length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips)]);
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
IsoD1MatrixVelX = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixVelY = NaN(minFramesInFlip,maxFlips,numFlies,2);
IsoD1MatrixVelTheta = NaN(minFramesInFlip,maxFlips,numFlies,2);
% 3D matrices with following structure: (flip #, fly#, odd/even flip)
IsoD1MatrixUpPCInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixUpGCInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixUpCircInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixUpRetInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixDownPCInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixDownGCInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixDownCircInFlip = NaN(maxFlips,numFlies,2);
IsoD1MatrixDownRetInFlip = NaN(maxFlips,numFlies,2);

% Now fill in the matrices we initialized above with all the data we have
for flyCounter = 1:numFlies
    % Do it for odd flips
    for oddFlipCounter = 1:length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips)
        % We're only plotting the data from frames with at least as many
        % frames as minFramesInFlip, so don't change the NaN if it doesn't
        % meet that criteria
        if length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).FrameInFlip) >= minFramesInFlip+1
            % Grab all the locomotion data
            IsoD1MatrixTime(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedTimeInFlip(1:minFramesInFlip);
            IsoD1MatrixX(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedCentroidX(1:minFramesInFlip);
            IsoD1MatrixY(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedCentroidY(1:minFramesInFlip);
            IsoD1MatrixTheta(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVerticalTheta(1:minFramesInFlip);
            IsoD1MatrixVelX(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelX(1:minFramesInFlip);
            IsoD1MatrixVelY(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelY(1:minFramesInFlip);
            IsoD1MatrixVelTheta(:,oddFlipCounter,flyCounter,1) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).AlignedFoldedVelTheta(1:minFramesInFlip);
            % Grab all the info about what crossing events happened in flip
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).UpProperCrosses) > 0
                IsoD1MatrixUpPCInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixUpPCInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).UpGlassCrosses) > 0
                IsoD1MatrixUpGCInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixUpGCInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).UpCircumventions) > 0
                IsoD1MatrixUpCircInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixUpCircInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).UpRetreats) > 0
                IsoD1MatrixUpRetInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixUpRetInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).DownProperCrosses) > 0
                IsoD1MatrixDownPCInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixDownPCInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).DownGlassCrosses) > 0
                IsoD1MatrixDownGCInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixDownGCInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).DownCircumventions) > 0
                IsoD1MatrixDownCircInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixDownCircInFlip(oddFlipCounter,flyCounter,1) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.OddFlips(oddFlipCounter).DownRetreats) > 0
                IsoD1MatrixDownRetInFlip(oddFlipCounter,flyCounter,1) = 1;
            else
                IsoD1MatrixDownRetInFlip(oddFlipCounter,flyCounter,1) = 0;
            end

        end
    end
    % Do it for even flips
    for evenFlipCounter = 1:length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips)
        % We're only plotting the data from frames with at least as many
        % frames as minFramesInFlip, so don't change the NaN if it doesn't
        % meet that criteria
        if length(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).FrameInFlip) >= minFramesInFlip+1
            IsoD1MatrixTime(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedTimeInFlip(1:minFramesInFlip);
            IsoD1MatrixX(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedCentroidX(1:minFramesInFlip);
            IsoD1MatrixY(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedCentroidY(1:minFramesInFlip);
            IsoD1MatrixTheta(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVerticalTheta(1:minFramesInFlip);
            IsoD1MatrixVelX(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelX(1:minFramesInFlip);
            IsoD1MatrixVelY(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelY(1:minFramesInFlip);
            IsoD1MatrixVelTheta(:,evenFlipCounter,flyCounter,2) = ...
                WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).AlignedFoldedVelTheta(1:minFramesInFlip);
            % Grab all the info about what crossing events happened in flip
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).UpProperCrosses) > 0
                IsoD1MatrixUpPCInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixUpPCInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).UpGlassCrosses) > 0
                IsoD1MatrixUpGCInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixUpGCInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).UpCircumventions) > 0
                IsoD1MatrixUpCircInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixUpCircInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).UpRetreats) > 0
                IsoD1MatrixUpRetInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixUpRetInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).DownProperCrosses) > 0
                IsoD1MatrixDownPCInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixDownPCInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).DownGlassCrosses) > 0
                IsoD1MatrixDownGCInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixDownGCInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).DownCircumventions) > 0
                IsoD1MatrixDownCircInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixDownCircInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
            if sum(WS.FlipBinnedFlyStruct(flyCounter).AlignedData.EvenFlips(evenFlipCounter).DownRetreats) > 0
                IsoD1MatrixDownRetInFlip(evenFlipCounter,flyCounter,2) = 1;
            else
                IsoD1MatrixDownRetInFlip(evenFlipCounter,flyCounter,2) = 0;
            end
        end
    end
end

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
plot(IsoD1MatrixTime(:,startsLowOddEven),...
         IsoD1MatrixY(:,startsLowOddEven),...
         'Color',[0.7,0.7,0.7])

% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
hold on
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
chosenFlyNum = 1;
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