% Pick which fly you want with flyNum, choose which odd flip you want to
% look at with oddFlipNum, and input what event you're looking to find
% using Elizabeth's eventCode scheme
% Run the top section to find the frame within the flip during which the
% fly initiates that gap event, then look for that frame number inside of
% UniqCompIdIndex to figure out what frameRange should be
% Set frameRange to be the range of frames you want to plot in the second
% section and just change the color being plotted depending on eventCode
% Section 3 allows you to overlay the skeleton on your plots

%% Section 1, finding the frames within flip
flyNum = 2;
oddFlipNum = 21;
eventCode = 440;
find(data(flyNum).BehavData.OddFlips(oddFlipNum).EventSeq == eventCode)

%% Section 2, plotting
figure
hold on
frameRange = 125:210;
% The 253 is approximately how many pixels separate the x coordinates of
% the centers of each corridor within a cassette
% The 185 is roughly the x coordinate in pixels of the center of the first
% corridor
% By doing these operations with 253 and 185, all flies and events should
% be able to be plotted on top of each other, regardless of corridor and
% independent of which side of the corridor the fly was on
% The 900 lets us flip the plots to match the real video (since the origin
% of the video is the top left corner)
if (data(flyNum).BehavData.OddFlips(oddFlipNum).CentroidX(frameRange(1)) - 253*(flyNum-1)) > 185
    plot((377-(data(flyNum).BehavData.OddFlips(oddFlipNum).CentroidX(frameRange) - 253*(flyNum-1))),900-data(flyNum).BehavData.OddFlips(oddFlipNum).CentroidY(frameRange),'b.');
else
    plot((data(flyNum).BehavData.OddFlips(oddFlipNum).CentroidX(frameRange) - 253*(flyNum-1)),900-data(flyNum).BehavData.OddFlips(oddFlipNum).CentroidY(frameRange),'b.');
end
xlim([75, 185]); % Approximate x coordinate limits of half the corridor
ylim([515, 620]); % Approximate y coordinate limits around 2.5 mm gap
title(oddFlipNum); % Helps track what you've plotted by giving a title

%% Section 3, overlaying the skeleton
hold on
figure
plot(WS.Skeleton_x([13, 13, 15, 15, 13, 13]),900-WS.Skeleton_y([13,14,14,17,17,18]), 'k--')
xlim([70, 185]);
ylim([515, 620]);