t1 = tiledlayout(1,7,"TileSpacing","tight");

% VisualizeRawDataScripts;
% close all

% Plot Y vs X in a histogram for all flies
for oddEven = 1:2
    nexttile()
    numBins = 400;
    histX = reshape([IsoD1MatrixX(:,:,:,oddEven)],[],1);
    histX = histX(~isnan(histX));
    histY = reshape([IsoD1MatrixY(:,:,:,oddEven)],[],1);
    histY = histY(~isnan(histY));
    XBinLim = 10;
    YBinLim = 25;
    histogram2(histX,histY,-XBinLim:(XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
    set(gca,'ColorScale','log');
    % colorbar;
    grid off
    % Example of how to overlay the skeleton in the plots
    xl = xlim;
    yl = ylim;
    hold on
    horizShift = 0.15;
    % a = plot([0;StraightenedSkelX(end/2+1:end,1);0],...
    %          [StraightenedSkelY(19,1);StraightenedSkelY(end/2+1:end,1);StraightenedSkelY(36,1)],'k');
    a = plot([StraightenedSkelX(:,oddEven);StraightenedSkelX(1,oddEven)]-horizShift,...
             [StraightenedSkelY(:,oddEven);StraightenedSkelY(1,oddEven)],'k');
    uistack(a,'top');
    uistack(a,'top');
    hold off
    ylim([-25,25]);
    % xlabel('X Position (mm)');
    % ylabel('Y Position (mm)');
    cmap = colormap('gray');
    colormap(flip(cmap,1));
    pbaspect([2*XBinLim 2*YBinLim 1])
    set(gca, 'color', 'none');
    set(gca,'visible','off')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

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
nexttile([1,3]);
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
     nanmedian(IsoD1MatrixY(:,startsLowOddEven),2),'k','LineWidth',3)
randSample = false(size(startsLow));
randSample(randi(size(startsLow,1),[1,5]),randi(size(startsLow,2),[1,5]),oddEven) = true;
randSample = and(randSample,startsLowOddEven);

prevRand = zeros(101*28*2,1);
prevRand([636,674,1343]) = 1;
prevRand = logical(reshape(prevRand,size(randSample)));

% Set this to 1 if you want to use previously randomly selected trajectories
prevVals = 1;

if prevVals == 1
    plot(IsoD1MatrixTime(:,prevRand),...
             IsoD1MatrixY(:,prevRand),...
             'b','LineWidth',1.5)
else
    plot(IsoD1MatrixTime(:,randSample),...
             IsoD1MatrixY(:,randSample),...
             'b','LineWidth',1.5)
end

% Example of how to overlay the gap locations in the plots
xl = xlim;
yl = ylim;
for i = [1,2,5,6,9,10,13,14,17,18]
    plot(xl(1):(xl(2)-xl(1))/10:xl(2),SkelY(i,oddEven)*ones(11,1),'k--')
end
hold off
xlim([0,max(IsoD1MatrixTime(:,randSample),[],'all')])

pbaspect([3.7*max(IsoD1MatrixTime(:,randSample),[],'all')/(2*XBinLim),1,1])
% title('y(t) for some flips starting below certain y position of all flies')
xlabel('Time in flip (sec)')
ylabel('Y position (mm)')
set(gca, 'color', 'none');
% set(gca,'visible','off')
% set(gca,'xtick',[])
% set(gca,'ytick',[])


% Plot one random good trajectory over time

flyNum = randi(size(IsoD1MatrixX,3)); % 2,19,15
flipNum = randi(size(IsoD1MatrixX,2)); % 18,56,3
flyNum = 2; % Selected this by sifting through random trajectories
flipNum = 18; % Selected this by sifting through random trajectories
oddEven = 1;


nexttile();
hold on
for i = 1:(size(IsoD1MatrixX,1)-1)
    plot(IsoD1MatrixTime(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixY(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
end
hold off
xlabel('Time (sec)');
ylabel('Y Position (mm)');
ylim([-25,25])
yt = gca;
pbaspect([20, 2*40, 1])
set(gca, 'color', 'none');
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])

XBinLim = 10;

nexttile();
hold on
for i = 1:(size(IsoD1MatrixX,1)-1)
    plot(IsoD1MatrixX(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixY(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
end
title(['Flip Num ', num2str(flipNum), ', Fly Num ', num2str(flyNum)])
pbaspect([2*XBinLim, 2*YBinLim, 1])

% a = plot([0;StraightenedSkelX(end/2+1:end,1);0],...
%          [StraightenedSkelY(19,1);StraightenedSkelY(end/2+1:end,1);StraightenedSkelY(36,1)],'k');
a = plot([StraightenedSkelX(:,oddEven);StraightenedSkelX(1,oddEven)]-horizShift,...
         [StraightenedSkelY(:,oddEven);StraightenedSkelY(1,oddEven)],'k');
uistack(a,'top');
uistack(a,'top');
hold off
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
traj = gca;
set(gca, 'color', 'none');
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([-XBinLim,XBinLim])
ylim([-25,25])
sizeOfFinalPanel = get(gca,'Position');
set(gcf,'Position',[4,384,1920,521]) % Chose these via trial and error for correct aspect ratios and sizes 


figure
set(gca, 'color', 'none');
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
hold on
for i = 1:(size(IsoD1MatrixX,1)-1)
    plot(IsoD1MatrixX(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixTime(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
end
hold off
xlabel('X Position (mm)');
ylabel('Time (sec)');
xlim([-XBinLim,XBinLim])
% pbaspect([2*XBinLim,IsoD1MatrixTime(size(IsoD1MatrixX,1),flipNum,flyNum,oddEven)/(2*YBinLim),1])
xt = gca;
set(gca,'Position',sizeOfFinalPanel)
set(gcf,'Position',[4,384,1920,521]) % Chose these via trial and error for correct aspect ratios and sizes 
% pbaspect([10, 8, 1])