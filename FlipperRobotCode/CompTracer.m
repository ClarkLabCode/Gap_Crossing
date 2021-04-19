% Asks user to outline the compartments and creates masks out of them

% Inputs
% inputFileName     = Name of video file to be analyzed, must include extension
% NumCorridors      = Number of corridors in the cassette
% NumGaps           = Number of gaps within each corridor
% indPos            = Index of frame at which flips happen

% Outputs
% CompMask1     = Cell array that holds the mask for each compartment for all
%                 odd numbered flips
% CompMask2     = Cell array that holds the mask for each compartment for all
%                 even numbered flips
% NumComps      = Number of compartments per corridor (either 3*NumGaps+1
%                 or 4*NumGaps+1, see below for more info)

% These use the compartment labeling scheme of 1 to 2*NumGaps+1 for the
% central parts of the corridors, then 2*NumGaps+2 to 3*NumGaps+1 for the
% well parts of the corridors

function [CompMask1, CompMask2, NumComps] = CompTracer(inputFileName, NumCorridors, NumGaps, indPos)

% Declare VideoReader object
reader1 = VideoReader(inputFileName);

% If NumComps = 3*NumGaps+1, then there is no distinction between left and
% right wells (the regions in which circumventions take place)
% If NumComps = 4*NumGaps+1, then there is a distinction between left/right
NumComps = 3*NumGaps+1;

% As in CorrMask, these cell arrays below hold masks and 1/2
% correspomnds to the parity of the flip number (1 = odd, 2 = even)

% CompMask is the final cell array that contains all masks of interest
CompMask1 = cell(1, NumComps);
CompMask2 = cell(1, NumComps);

% tempCentCorrMask is the cell array that is used as an intermediate to
% label the central parts of each corridor
tempCentCorrMask1 = cell(1, NumCorridors);
tempCentCorrMask2 = cell(1, NumCorridors);

% tempGapMask is the cell array that is used as an intermediate to label
% the gap parts of each corridor
tempGapMask1 = cell(1, NumGaps);
tempGapMask2 = cell(1, NumGaps);

% Let the user draw the mask for the center of each corridor and save
% within tempCentCorrMask1/2. Gives text output to ensure user keeps track
% of which corridor they must label.
for i = 1:NumCorridors
    fprintf(strcat('Select Center Box of Corridor #',num2str(i),' in orientation 1\n'));
    tempCentCorrMask1{i} = roipoly(read(reader1,indPos(1)+5));
    fprintf(strcat('Select Center Box of Corridor #',num2str(i),' in orientation 2\n'));
    tempCentCorrMask2{i} = roipoly(read(reader1,indPos(3)+5)); 
end

% Let the user draw the mask for all the gaps of the same width and save
% within tempGapMask1/2. Gives text output to ensure user keeps track
% of which gap size they must label.
for i = 1:NumGaps
    fprintf('Select Gap Box %.1fmm in orientation 1\n',(2+0.5*i));
    tempGapMask1{i} = roipoly(read(reader1,indPos(1)+5));
    fprintf('Select Gap Box %.1fmm in orientation 2\n',(2+0.5*i));
    tempGapMask2{i} = roipoly(read(reader1,indPos(3)+5)); 
end

% Defining some other cell arrays that are used to compute the comps
% Doing it this way ensures the full tiling of space in the video and also
% minimizes the amount of user inputs required substantially
tempAntiGapMask1 = cell(1, NumGaps+1);
tempAntiGapMask2 = cell(1, NumGaps+1);
AntiGapMask1 = cell(1, NumGaps+1);
AntiGapMask2 = cell(1, NumGaps+1);
GapMask1 = cell(1, NumGaps);
GapMask2 = cell(1, NumGaps);
tempAntiMask1 = true(size(tempGapMask1{1},1),size(tempGapMask1{1},2));
tempAntiMask2 = true(size(tempGapMask2{1},1),size(tempGapMask1{1},2));
tempVertCentMask1 = false(size(tempGapMask1{1},1),size(tempGapMask1{1},2));
tempVertCentMask2 = false(size(tempGapMask2{1},1),size(tempGapMask1{1},2));

% tempAnitMask is a single logical array that is 1 everywhere except in the
% regions that the user identified gaps
for i = 1:NumGaps
    tempAntiMask1 = logical(tempAntiMask1-tempGapMask1{i});
    tempAntiMask2 = logical(tempAntiMask2-tempGapMask2{i});
end

% tempVertCentMask is a single logical array that is 1 only in the regions that
% the user identified as being the centers of each corridor
for i = 1:NumCorridors
    tempVertCentMask1 = logical(tempVertCentMask1+tempCentCorrMask1{i});
    tempVertCentMask2 = logical(tempVertCentMask1+tempCentCorrMask2{i});
end

% NonGapMask is a single logical array that is 1 only in the regions that
% the user identified as being in the center of each corridor but not
% within a gap
NonGapMask1 = logical(tempVertCentMask1.*tempAntiMask1);
NonGapMask2 = logical(tempVertCentMask2.*tempAntiMask2);

% Let the user draw the mask for all the regions in which there are no gaps
% and save within tempAntiGapMask1/2. Gives text output to ensure user 
% keeps track of which non-gap size they must label. 
for i = 1:(NumGaps+1)
    fprintf('Select Non-Gap Box %.1fmm in orientation 1.\nStop roughly at the gap.\n',(2+0.5*i));
    tempAntiGapMask1{i} = roipoly(read(reader1,indPos(1)+5));
    fprintf('Select Non-Gap Box %.1fmm in orientation 2.\nStop roughly at the gap.\n',(2+0.5*i));
    tempAntiGapMask2{i} = roipoly(read(reader1,indPos(3)+5)); 
end

% AntiGapMask is a cell array that contains the mask corresponding to each
% non-gap region of the corridors
% e.g. AntiGapMask{1} is the mask corresponding to the region of all the
% corridors that is above or below (depending on orientation) the 2.5mm gap
% and AntiGapMask{2} is the mask corresponding to the region of all the
% corridors that is between the 2.5mm and 3mm gaps
for i = 1:(NumGaps+1)
    AntiGapMask1{i} = logical(tempAntiGapMask1{i}.*NonGapMask1);
    AntiGapMask2{i} = logical(tempAntiGapMask2{i}.*NonGapMask2);
end

% Gap Mask is a cell array that contains the mask corresponding to each gap
% region of the corridors, excluding the wells
% e.g. GapMask{1} is the mask corresponding to the region of all the
% corridors that contrains the 2.5mm gap but not the wells
for i = 1:NumGaps
    GapMask1{i} = logical(tempVertCentMask1.*tempGapMask1{i});
    GapMask2{i} = logical(tempVertCentMask2.*tempGapMask2{i});
end

% Define cell array that will hold the well masks
WellMask1 = cell(1, NumGaps);
WellMask2 = cell(1, NumGaps);

% WellMask is a cell array that contains the mask corresponding to only the
% wells within each gap region of the corridors
% e.g. WellMask{1} is the mask corresponding to the region of all the
% corridors that contains only the wells of the 2.5mm gap
for i = 1:NumGaps
    WellMask1{i} = logical(tempGapMask1{i} - GapMask1{i});
    WellMask2{i} = logical(tempGapMask2{i} - GapMask2{i});
end

% Now that we have all the necessary masks, fill up CompMask1/2
% This uses the compartment labeling scheme of 1 to 2*NumGaps+1 for the
% central parts of the corridors, then 2*NumGaps+2 to 3*NumGaps+1 for the
% well parts of the corridors
for i = 1:(NumGaps+1)
    CompMask1{2*i-1} = AntiGapMask1{i};
    CompMask2{2*i-1} = AntiGapMask2{i};
end

for i = 1:NumGaps
    CompMask1{2*i} = GapMask1{i};
    CompMask2{2*i} = GapMask2{i};
    CompMask1{i+2*NumGaps+1} = WellMask1{i};
    CompMask2{i+2*NumGaps+1} = WellMask2{i};
end

end
