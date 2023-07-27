% Computes and plots all gap crossing-related statistics
% Breaks down for even/odd flips and for up/down motion

% "Crossing" Rate vs Gap Width
% Frac of Crossings that are Proper Crossings vs Gap Width
% Proper Crossing Rate vs Gap Width

function WS = AnalyzeCrossStats(WS)

% Make sure there are no conflicting figures open
close all

% Move to the correct data folder so that files can be saved in right place
flipperRobotCodePath = pwd;
dataPath = [flipperRobotCodePath,'\..\..\Data\'];
cd(dataPath);
dataPath = pwd;
dataPath = [dataPath, '\'];
cd(flipperRobotCodePath);
% If the data being analyzed came from combining mismatched experiments,
% ask the user where to save the file to/from
if strcmp(WS.directoryName,'Potentially variable')
    warning('Beware, analyzing data that came from combining mismatched experiments.')
    directoryName = uigetdir(dataPath,'Select folder to save this analyzed WS');
    cd(directoryName);
else
    cd([dataPath, WS.directoryName]);
end

% Port in the relevant fields from WS
FBFS = WS.FlipBinnedFlyStruct;
NumGaps = WS.NumGaps;

% This gives this function the ability to be used on single experiments or
% on structs that contain multiple experiments since the FBFS struct
% contains a dummy ExpNum layer when it is only one experiment but that
% layer doesn't exist when there are multiple experiments combined by
% CombineExpStructs
if isfield(FBFS,'ExpNum')
    % Just a sanity check that there isn't something very weird with FBFS
    if length(FBFS) ~= 1
        error('Something is wrong with your FlipBinnedFlyStruct')
    end
    % Remove the dummy ExpNum layer
    FBFS = FBFS.ExpNum;
end

% Find number of flies in the struct
NumFlies = length(FBFS);

% Initialize matrices that hold # of gap events for each width, for each fly
% Have such a matrix for every combination of up/down and even/odd

% Up + Odd
UpOddVecCircumvent      = zeros(NumFlies,NumGaps);
UpOddVecRetreat         = zeros(NumFlies,NumGaps);
UpOddVecGlassCross      = zeros(NumFlies,NumGaps);
UpOddVecProperCross     = zeros(NumFlies,NumGaps);

% Down + Odd
DownOddVecCircumvent    = zeros(NumFlies,NumGaps);
DownOddVecRetreat       = zeros(NumFlies,NumGaps);
DownOddVecGlassCross    = zeros(NumFlies,NumGaps);
DownOddVecProperCross   = zeros(NumFlies,NumGaps);

% Up + Even
UpEvenVecCircumvent     = zeros(NumFlies,NumGaps);
UpEvenVecRetreat        = zeros(NumFlies,NumGaps);
UpEvenVecGlassCross     = zeros(NumFlies,NumGaps);
UpEvenVecProperCross    = zeros(NumFlies,NumGaps);

% Down + Even
DownEvenVecCircumvent   = zeros(NumFlies,NumGaps);
DownEvenVecRetreat      = zeros(NumFlies,NumGaps);
DownEvenVecGlassCross   = zeros(NumFlies,NumGaps);
DownEvenVecProperCross  = zeros(NumFlies,NumGaps);

% Go through all the flies and fill up these matrices
% Note that we're summing across all flips per fly at each width
% (this can be seen by the fact that the matrix element appears on both the
% right and left hand side of the equality when looping through flips)
for flyCounter = 1:NumFlies
    for gapCounter = 1:NumGaps
        % Odd flips
        for oddFlipCounter = 1:length(FBFS(flyCounter).BehavData.OddFlips)
            % Up + Odd
            UpOddVecCircumvent(flyCounter,gapCounter) = ...
                UpOddVecCircumvent(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).UpCircumventions(gapCounter);
            UpOddVecRetreat(flyCounter,gapCounter) = ...
                UpOddVecRetreat(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).UpRetreats(gapCounter);
            UpOddVecGlassCross(flyCounter,gapCounter) = ...
                UpOddVecGlassCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).UpGlassCrosses(gapCounter);
            UpOddVecProperCross(flyCounter,gapCounter) = ...
                UpOddVecProperCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).UpProperCrosses(gapCounter);
            % Down + Odd
            DownOddVecCircumvent(flyCounter,gapCounter) = ...
                DownOddVecCircumvent(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).DownCircumventions(gapCounter);
            DownOddVecRetreat(flyCounter,gapCounter) = ...
                DownOddVecRetreat(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).DownRetreats(gapCounter);
            DownOddVecGlassCross(flyCounter,gapCounter) = ...
                DownOddVecGlassCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).DownGlassCrosses(gapCounter);
            DownOddVecProperCross(flyCounter,gapCounter) = ...
                DownOddVecProperCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.OddFlips(oddFlipCounter).DownProperCrosses(gapCounter);
        end
        
        % Even flips
        for evenFlipCounter = 1:length(FBFS(flyCounter).BehavData.EvenFlips)
            % Up + Even
            UpEvenVecCircumvent(flyCounter,gapCounter) = ...
                UpEvenVecCircumvent(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpCircumventions(gapCounter);
            UpEvenVecRetreat(flyCounter,gapCounter) = ...
                UpEvenVecRetreat(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpRetreats(gapCounter);
            UpEvenVecGlassCross(flyCounter,gapCounter) = ...
                UpEvenVecGlassCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpGlassCrosses(gapCounter);
            UpEvenVecProperCross(flyCounter,gapCounter) = ...
                UpEvenVecProperCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).UpProperCrosses(gapCounter);
            % Down + Even
            DownEvenVecCircumvent(flyCounter,gapCounter) = ...
                DownEvenVecCircumvent(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownCircumventions(gapCounter);
            DownEvenVecRetreat(flyCounter,gapCounter) = ...
                DownEvenVecRetreat(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownRetreats(gapCounter);
            DownEvenVecGlassCross(flyCounter,gapCounter) = ...
                DownEvenVecGlassCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownGlassCrosses(gapCounter);
            DownEvenVecProperCross(flyCounter,gapCounter) = ...
                DownEvenVecProperCross(flyCounter,gapCounter) + ...
                FBFS(flyCounter).BehavData.EvenFlips(evenFlipCounter).DownProperCrosses(gapCounter);
        end
    end
end

% Now start making vectors that hold information about just even/odd or
% just up/down by summing over the opposite pair

% Look at all odd by summing (odd + up) and (odd + down)
AllOddVecCircumvent     = [UpOddVecCircumvent    + DownOddVecCircumvent];
AllOddVecRetreat        = [UpOddVecRetreat       + DownOddVecRetreat];
AllOddVecGlassCross     = [UpOddVecGlassCross    + DownOddVecGlassCross];
AllOddVecProperCross    = [UpOddVecProperCross   + DownOddVecProperCross];

AllOddVecAllGapEvents        = ...
    AllOddVecCircumvent + AllOddVecRetreat + AllOddVecGlassCross + AllOddVecProperCross;
AllOddVecAllCrossEvents      = ...
    AllOddVecGlassCross + AllOddVecProperCross;

% Look at all even by summing (even + up) and (even + down)
AllEvenVecCircumvent    = [UpEvenVecCircumvent   + DownEvenVecCircumvent];
AllEvenVecRetreat       = [UpEvenVecRetreat      + DownEvenVecRetreat];
AllEvenVecGlassCross    = [UpEvenVecGlassCross   + DownEvenVecGlassCross];
AllEvenVecProperCross   = [UpEvenVecProperCross  + DownEvenVecProperCross];

AllEvenVecAllGapEvents        = ...
    AllEvenVecCircumvent + AllEvenVecRetreat + AllEvenVecGlassCross + AllEvenVecProperCross;
AllEvenVecAllCrossEvents      = ...
    AllEvenVecGlassCross + AllEvenVecProperCross;

% Look at all up by summing (odd + up) and (even + up)
AllUpVecCircumvent      = [UpOddVecCircumvent    + UpEvenVecCircumvent];
AllUpVecRetreat         = [UpOddVecRetreat       + UpEvenVecRetreat];
AllUpVecGlassCross      = [UpOddVecGlassCross    + UpEvenVecGlassCross];
AllUpVecProperCross     = [UpOddVecProperCross   + UpEvenVecProperCross];

AllUpVecAllGapEvents        = ...
    AllUpVecCircumvent + AllUpVecRetreat + AllUpVecGlassCross + AllUpVecProperCross;
AllUpVecAllCrossEvents      = ...
    AllUpVecGlassCross + AllUpVecProperCross;

% Look at all down by summing (odd + down) and (even + down)
AllDownVecCircumvent    = [DownOddVecCircumvent  + DownEvenVecCircumvent];
AllDownVecRetreat       = [DownOddVecRetreat     + DownEvenVecRetreat];
AllDownVecGlassCross    = [DownOddVecGlassCross  + DownEvenVecGlassCross];
AllDownVecProperCross   = [DownOddVecProperCross + DownEvenVecProperCross];

AllDownVecAllGapEvents        = ...
    AllDownVecCircumvent + AllDownVecRetreat + AllDownVecGlassCross + AllDownVecProperCross;
AllDownVecAllCrossEvents      = ...
    AllDownVecGlassCross + AllDownVecProperCross;

% Look at all by summing over (all up) and (all down)
AllAllVecCircumvent         = [AllUpVecCircumvent    + AllDownVecCircumvent];
AllAllVecRetreat            = [AllUpVecRetreat       + AllDownVecRetreat];
AllAllVecGlassCross         = [AllUpVecGlassCross    + AllDownVecGlassCross];
AllAllVecProperCross        = [AllUpVecProperCross   + AllDownVecProperCross];

AllAllVecAllGapEvents = ...
    AllAllVecCircumvent + AllAllVecRetreat + AllAllVecGlassCross + AllAllVecProperCross;
AllAllVecAllCrossEvents = ...
    AllAllVecGlassCross + AllAllVecProperCross;

% Now generate the rate vectors that we're interested in
% 1st line: "Crossings" / All Gap Events
% 2nd line: Proper Crossings / "Crossings"
% 3rd line: Proper Crossings / All Gap Events

% Do it for all odd
AllOddVecAllCrossOverAllGapEventRate        = AllOddVecAllCrossEvents./AllOddVecAllGapEvents;
AllOddVecProperCrossOverAllCrossRate        = AllOddVecProperCross./AllOddVecAllCrossEvents;
AllOddVecProperCrossOverAllGapEventRate     = AllOddVecProperCross./AllOddVecAllGapEvents;

% Do it for all even
AllEvenVecAllCrossOverAllGapEventRate       = AllEvenVecAllCrossEvents./AllEvenVecAllGapEvents;
AllEvenVecProperCrossOverAllCrossRate       = AllEvenVecProperCross./AllEvenVecAllCrossEvents;
AllEvenVecProperCrossOverAllGapEventRate    = AllEvenVecProperCross./AllEvenVecAllGapEvents;

% Do it for all up
AllUpVecAllCrossOverAllGapEventRate         = AllUpVecAllCrossEvents./AllUpVecAllGapEvents;
AllUpVecProperCrossOverAllCrossRate         = AllUpVecProperCross./AllUpVecAllCrossEvents;
AllUpVecProperCrossOverAllGapEventRate      = AllUpVecProperCross./AllUpVecAllGapEvents;

% Do it for all down
AllDownVecAllCrossOverAllGapEventRate       = AllDownVecAllCrossEvents./AllDownVecAllGapEvents;
AllDownVecProperCrossOverAllCrossRate       = AllDownVecProperCross./AllDownVecAllCrossEvents;
AllDownVecProperCrossOverAllGapEventRate    = AllDownVecProperCross./AllDownVecAllGapEvents;

% Do it for all (i.e., don't care if up/down and even/odd)
AllAllVecAllCrossOverAllGapEventRate        = AllAllVecAllCrossEvents./AllAllVecAllGapEvents;
AllAllVecProperCrossOverAllCrossRate        = AllAllVecProperCross./AllAllVecAllCrossEvents;
AllAllVecProperCrossOverAllGapEventRate     = AllAllVecProperCross./AllAllVecAllGapEvents;

% The sizes of the gaps in the experiments (in mm)
gapSizes = 1:0.5:2.5;

% Now we want to filter out any flies that didn't perform at least 10 gap
% events at every width
% Similarly, filter out flies that didn't perform at least 10 "crossings" 
% at all gap widths for the calculation of Proper Crossings / "Crossings"

% We do this by taking the minimum number of events across widths and
% checking if it is at least 10 (hence the min( __ , [],2) >= 10) and then
% using this logical vector to index the flies that do have at least 10 events

% All odd
% Generate logical vectors
AllOddVecAllGapEventsLog = min(AllOddVecAllGapEvents,[],2)>=10;
AllOddVecAllCrossEventsLog = min(AllOddVecAllCrossEvents,[],2)>=10;

% Use the logical vectors to index the full data
AllOddVecAllCrossOverAllGapEventRate = ...
    AllOddVecAllCrossOverAllGapEventRate(AllOddVecAllGapEventsLog,:);
AllOddVecProperCrossOverAllCrossRate = ...
    AllOddVecProperCrossOverAllCrossRate(AllOddVecAllCrossEventsLog,:);
AllOddVecProperCrossOverAllGapEventRate = ...
    AllOddVecProperCrossOverAllGapEventRate(AllOddVecAllGapEventsLog,:);

% All even
% Generate logical vectors
AllEvenVecAllGapEventsLog = min(AllEvenVecAllGapEvents,[],2)>=10;
AllEvenVecAllCrossEventsLog = min(AllEvenVecAllCrossEvents,[],2)>=10;

% Use the logical vectors to index the full data
AllEvenVecAllCrossOverAllGapEventRate = ...
    AllEvenVecAllCrossOverAllGapEventRate(AllEvenVecAllGapEventsLog,:);
AllEvenVecProperCrossOverAllCrossRate = ...
    AllEvenVecProperCrossOverAllCrossRate(AllEvenVecAllCrossEventsLog,:);
AllEvenVecProperCrossOverAllGapEventRate = ...
    AllEvenVecProperCrossOverAllGapEventRate(AllEvenVecAllGapEventsLog,:);

% All up
% Generate logical vectors
AllUpVecAllGapEventsLog = min(AllUpVecAllGapEvents,[],2)>=10;
AllUpVecAllCrossEventsLog = min(AllUpVecAllCrossEvents,[],2)>=10;

% Use the logical vectors to index the full data
AllUpVecAllCrossOverAllGapEventRate = ...
    AllUpVecAllCrossOverAllGapEventRate(AllUpVecAllGapEventsLog,:);
AllUpVecProperCrossOverAllCrossRate = ...
    AllUpVecProperCrossOverAllCrossRate(AllUpVecAllCrossEventsLog,:);
AllUpVecProperCrossOverAllGapEventRate = ...
    AllUpVecProperCrossOverAllGapEventRate(AllUpVecAllGapEventsLog,:);

% All down
% Generate logical vectors
AllDownVecAllGapEventsLog = min(AllDownVecAllGapEvents,[],2)>=10;
AllDownVecAllCrossEventsLog = min(AllDownVecAllCrossEvents,[],2)>=10;

% Use the logical vectors to index the full data
AllDownVecAllCrossOverAllGapEventRate = ...
    AllDownVecAllCrossOverAllGapEventRate(AllDownVecAllGapEventsLog,:);
AllDownVecProperCrossOverAllCrossRate = ...
    AllDownVecProperCrossOverAllCrossRate(AllDownVecAllCrossEventsLog,:);
AllDownVecProperCrossOverAllGapEventRate = ...
    AllDownVecProperCrossOverAllGapEventRate(AllDownVecAllGapEventsLog,:);

% All events (i.e., don't care if up/down and even/odd)
% Generate logical vectors
AllAllVecAllGapEventsLog = min(AllAllVecAllGapEvents,[],2)>=10;
AllAllVecAllCrossEventsLog = min(AllAllVecAllCrossEvents,[],2)>=10;

% Use the logical vectors to index the full data
AllAllVecAllCrossOverAllGapEventRate = ...
    AllAllVecAllCrossOverAllGapEventRate(AllAllVecAllGapEventsLog,:);
AllAllVecProperCrossOverAllCrossRate = ...
    AllAllVecProperCrossOverAllCrossRate(AllAllVecAllCrossEventsLog,:);
AllAllVecProperCrossOverAllGapEventRate = ...
    AllAllVecProperCrossOverAllGapEventRate(AllAllVecAllGapEventsLog,:);

% Now we compute the standard error for these rates computed above
% 1st line: "Crossings" / All Gap Events
% 2nd line: Proper Crossings / "Crossings"
% 3rd line: Proper Crossings / All Gap Events

% Do it for all events (i.e., don't care if up/down and even/odd)
AllAllVecAllCrossOverAllGapEventStd = ...
    std(AllAllVecAllCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllAllVecAllCrossOverAllGapEventRate));
AllAllVecProperCrossOverAllCrossStd = ...
    std(AllAllVecProperCrossOverAllCrossRate,0,1)...
        /sqrt(length(AllAllVecProperCrossOverAllCrossRate));
AllAllVecProperCrossOverAllGapEventStd = ...
    std(AllAllVecProperCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllAllVecProperCrossOverAllGapEventRate));
    
% All odd    
AllOddVecAllCrossOverAllGapEventStd = ...
    std(AllOddVecAllCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllOddVecAllCrossOverAllGapEventRate));
AllOddVecProperCrossOverAllCrossStd = ...
    std(AllOddVecProperCrossOverAllCrossRate,0,1)...
        /sqrt(length(AllOddVecProperCrossOverAllCrossRate));
AllOddVecProperCrossOverAllGapEventStd = ...
    std(AllOddVecProperCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllOddVecProperCrossOverAllGapEventRate));
    
% All even
AllEvenVecAllCrossOverAllGapEventStd = ...
    std(AllEvenVecAllCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllEvenVecAllCrossOverAllGapEventRate));
AllEvenVecProperCrossOverAllCrossStd = ...
    std(AllEvenVecProperCrossOverAllCrossRate,0,1)...
        /sqrt(length(AllEvenVecProperCrossOverAllCrossRate));
AllEvenVecProperCrossOverAllGapEventStd = ...
    std(AllEvenVecProperCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllEvenVecProperCrossOverAllGapEventRate));
    
% All up
AllUpVecAllCrossOverAllGapEventStd = ...
    std(AllUpVecAllCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllUpVecAllCrossOverAllGapEventRate));
AllUpVecProperCrossOverAllCrossStd = ...
    std(AllUpVecProperCrossOverAllCrossRate,0,1)...
        /sqrt(length(AllUpVecProperCrossOverAllCrossRate));
AllUpVecProperCrossOverAllGapEventStd = ...
    std(AllUpVecProperCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllUpVecProperCrossOverAllGapEventRate));
   
% All down
AllDownVecAllCrossOverAllGapEventStd = ...
    std(AllDownVecAllCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllDownVecAllCrossOverAllGapEventRate));
AllDownVecProperCrossOverAllCrossStd = ...
    std(AllDownVecProperCrossOverAllCrossRate,0,1)...
        /sqrt(length(AllDownVecProperCrossOverAllCrossRate));
AllDownVecProperCrossOverAllGapEventStd = ...
    std(AllDownVecProperCrossOverAllGapEventRate,0,1)...
        /sqrt(length(AllDownVecProperCrossOverAllGapEventRate));
    
% Now we begin to plot these curves
    
% figure
% 
% % Plot "Crossings" / All Events for even, odd, and all 
% hold on
% errorbar(gapSizes,mean(AllOddVecAllCrossOverAllGapEventRate,1),AllOddVecAllCrossOverAllGapEventStd);
% errorbar(gapSizes,mean(AllEvenVecAllCrossOverAllGapEventRate,1),AllEvenVecAllCrossOverAllGapEventStd);
% errorbar(gapSizes,mean(AllAllVecAllCrossOverAllGapEventRate,1),AllAllVecAllCrossOverAllGapEventStd);
% hold off
% 
% title([WS.directoryName, ' (', WS.genotype, ')'], 'Interpreter', 'none')
% xlabel('Gap Width (mm)');
% ylabel('All Cross Events / All Gap Events');
% xlim([0.8,2.7]);
% xticks(1:0.5:2.5);
% ylim([0,1]);
% legend('Odd Flips','Even Flips','All Flips');
%     
% figure
% 
% % Plot Proper Crossings / "Crossings" for even, odd, and all 
% hold on
% errorbar(gapSizes,mean(AllOddVecProperCrossOverAllCrossRate,1),AllOddVecProperCrossOverAllCrossStd);
% errorbar(gapSizes,mean(AllEvenVecProperCrossOverAllCrossRate,1),AllEvenVecProperCrossOverAllCrossStd);
% errorbar(gapSizes,mean(AllAllVecProperCrossOverAllCrossRate,1),AllAllVecProperCrossOverAllCrossStd);
% hold off
% 
% title([WS.directoryName, ' (', WS.genotype, ')'], 'Interpreter', 'none')
% xlabel('Gap Width (mm)');
% ylabel('Proper Cross Events / All Cross Events');
% xlim([0.8,2.7]);
% xticks(1:0.5:2.5);
% ylim([0,1]);
% legend('Odd Flips','Even Flips','All Flips');
%     
% figure
% 
% % Plot Proper Crossings / All Events for even, odd, and all 
% hold on
% errorbar(gapSizes,mean(AllOddVecProperCrossOverAllGapEventRate,1),AllOddVecProperCrossOverAllGapEventStd);
% errorbar(gapSizes,mean(AllEvenVecProperCrossOverAllGapEventRate,1),AllEvenVecProperCrossOverAllGapEventStd);
% errorbar(gapSizes,mean(AllAllVecProperCrossOverAllGapEventRate,1),AllAllVecProperCrossOverAllGapEventStd);
% hold off
% 
% title([WS.directoryName, ' (', WS.genotype, ')'], 'Interpreter', 'none')
% xlabel('Gap Width (mm)');
% ylabel('Proper Cross Events / All Gap Events');
% xlim([0.8,2.7]);
% xticks(1:0.5:2.5);
% ylim([0,1]);
% legend('Odd Flips','Even Flips','All Flips');

% figure
% 
% % Plot "Crossings" / All Events for up, down, and all 
% hold on
% errorbar(gapSizes,mean(AllUpVecAllCrossOverAllGapEventRate,1),AllUpVecAllCrossOverAllGapEventStd);
% errorbar(gapSizes,mean(AllDownVecAllCrossOverAllGapEventRate,1),AllDownVecAllCrossOverAllGapEventStd);
% errorbar(gapSizes,mean(AllAllVecAllCrossOverAllGapEventRate,1),AllAllVecAllCrossOverAllGapEventStd);
% hold off
% 
% title([WS.directoryName, ' (', WS.genotype, ')'], 'Interpreter', 'none')
% xlabel('Gap Width (mm)');
% ylabel('All Cross Events / All Gap Events');
% xlim([0.8,2.7]);
% xticks(1:0.5:2.5);
% ylim([0,1]);
% legend('Up','Down','All');
    
% figure
% 
% % Plot Proper Crossings / "Crossings" for up, down, and all 
% hold on
% errorbar(gapSizes,mean(AllUpVecProperCrossOverAllCrossRate,1),AllUpVecProperCrossOverAllCrossStd);
% errorbar(gapSizes,mean(AllDownVecProperCrossOverAllCrossRate,1),AllDownVecProperCrossOverAllCrossStd);
% errorbar(gapSizes,mean(AllAllVecProperCrossOverAllCrossRate,1),AllAllVecProperCrossOverAllCrossStd);
% hold off
% 
% title([WS.directoryName, ' (', WS.genotype, ')'], 'Interpreter', 'none')
% xlabel('Gap Width (mm)');
% ylabel('Proper Cross Events / All Cross Events');
% xlim([0.8,2.7]);
% xticks(1:0.5:2.5);
% ylim([0,1]);
% legend('Up','Down','All');
    
figure

% Plot Proper Crossings / All Events for up, down, and all 
hold on
errorbar(gapSizes,mean(AllUpVecProperCrossOverAllGapEventRate,1),AllUpVecProperCrossOverAllGapEventStd);
errorbar(gapSizes,mean(AllDownVecProperCrossOverAllGapEventRate,1),AllDownVecProperCrossOverAllGapEventStd);
errorbar(gapSizes,mean(AllAllVecProperCrossOverAllGapEventRate,1),AllAllVecProperCrossOverAllGapEventStd);
hold off

title([WS.directoryName, ' (', WS.genotype, ')'], 'Interpreter', 'none')
xlabel('Gap Width (mm)');
ylabel('Proper Cross Events / All Gap Events');
xlim([0.8,2.7]);
xticks(1:0.5:2.5);
ylim([0,1]);
legend('Up','Down','All');
       
% Clear the variables we don't care to save to WS
clear FBFS NumGaps NumFlies

% Populate a struct called crossStats that contains all the above computed
% quantities by first saving the workspace, loading it in as the struct
% crossStats, then remove crossStats.WS (since WS is part of the workspace)
% Finally, delete the dummy workspace file
save('dummyFileToBeDeleted.mat')
crossStats = load('dummyFileToBeDeleted.mat');
delete('dummyFileToBeDeleted.mat');
crossStats = rmfield(crossStats,'WS');

% Save the figures
% saveas(figure(1), [WS.directoryName, '_Crossing_Plot.fig']);
% saveas(figure(2), [WS.directoryName, '_Proper_Crossing_Plot.fig']);
saveas(figure(1), [WS.directoryName, '_Proper_Crossing_Plot.fig']);

% Update the fields in WS
WS.crossStats = crossStats;

% Now save the updated WS which has crossStats in it
save([WS.directoryName,'_WS_Combined.mat'],'WS')

% Navigate back to the data folder
cd(flipperRobotCodePath);

end
