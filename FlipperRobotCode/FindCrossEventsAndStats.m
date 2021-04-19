% Finds all crossings, retreats, and circumventions per fly and
% then plots the curves with error bars

% Beware, this function is not currently well-generalizable to more gaps

% Outputs:
% finalFlyStruct    = Fly-centric structure that holds all fly info
% meanCrossRate     = Gap crossing rate averaged across all flies
% stderror          = SEM for crossing rate (uses n = # flies)
% FlyCrossCountRate = Gap crossing rate for each individual fly
% FlyCrossBinomErr  = SEM for crossing rate for each fly (n = # of crossing events)

function [finalFlyStruct, meanCrossRate, stderror, FlyCrossCountRate, FlyCrossBinomErr, fig1, fig2] = ...
    FindCrossEventsAndStats(finalFlyStruct, NumGaps)

% Define the TransitionIDs for each gap crossing
Gap1CrossTransIDFor = [1,2,3];
Gap1CrossTransIDBack = [3,2,1];
Gap2CrossTransIDFor = [3,4,5];
Gap2CrossTransIDBack = [5,4,3];
Gap3CrossTransIDFor = [5,6,7];
Gap3CrossTransIDBack = [7,6,5];
Gap4CrossTransIDFor = [7,8,9];
Gap4CrossTransIDBack = [9,8,7];

% Define the TransitionIDs for each gap retreat
Gap1RetTransIDFor = [1,2,1];
Gap1RetTransIDBack = [3,2,3];
Gap2RetTransIDFor = [3,4,3];
Gap2RetTransIDBack = [5,4,5];
Gap3RetTransIDFor = [5,6,5];
Gap3RetTransIDBack = [7,6,7];
Gap4RetTransIDFor = [7,8,7];
Gap4RetTransIDBack = [9,8,9];

% Define the TransitionIDs for each gap circumvention
Gap1CircTransID = [2,10,2];
Gap2CircTransID = [4,11,4];
Gap3CircTransID = [6,12,6];
Gap4CircTransID = [8,13,8];

% Create vectors to hold cross, ret, circ info for each fly per gap width
FlyCrossCount = zeros(length(finalFlyStruct),NumGaps);
FlyRetCount = zeros(length(finalFlyStruct),NumGaps);
FlyCircCount = zeros(length(finalFlyStruct),NumGaps);

for flyStructCounter = 1:length(finalFlyStruct)
    % Find the frames in which a fly changes compartments
    CompTransVec = finalFlyStruct(flyStructCounter).CompID - ...
        circshift(finalFlyStruct(flyStructCounter).CompID,1);
    CompTransVec(1) = 0;
    
    % Find the frames in which the cassette had been flipped
    FlipTransVec = finalFlyStruct(flyStructCounter).FlipNumber - ...
        circshift(finalFlyStruct(flyStructCounter).FlipNumber,1);
    FlipTransVec(1) = 0;
   
    % Turn these into logical arrays
    CompTransLog = CompTransVec ~= 0;
    FlipTransLog = FlipTransVec ~= 0;
    
    % Don't count compartment changes between flips
    TrueCompTransLog = logical(CompTransLog.*(1-FlipTransLog));
    
    % Do I use this anymore?
%     GapEventIDs = CompTransVec(TrueCompTransLog);

    % Grab all the unique CompIDs from finalFlyStruct and index them
    CompTransIDVec = finalFlyStruct(flyStructCounter).CompID(TrueCompTransLog);
    finalFlyStruct(flyStructCounter).UniqCompID = CompTransIDVec;
    finalFlyStruct(flyStructCounter).UniqCompIDIndex = find(TrueCompTransLog);

    % Count number of times the TransIDs happen for gap crossings per width
    % and save the positions of the crossing events
    Gap1CrossCount = length(strfind(CompTransIDVec, Gap1CrossTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap1CrossTransIDBack));
    finalFlyStruct(flyStructCounter).Gap1CrossIndex = ...
        sort([strfind(CompTransIDVec, Gap1CrossTransIDFor) strfind(CompTransIDVec, Gap1CrossTransIDBack)]);
    Gap2CrossCount = length(strfind(CompTransIDVec, Gap2CrossTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap2CrossTransIDBack));
    finalFlyStruct(flyStructCounter).Gap2CrossIndex = ...
        sort([strfind(CompTransIDVec, Gap2CrossTransIDFor) strfind(CompTransIDVec, Gap2CrossTransIDBack)]);
    Gap3CrossCount = length(strfind(CompTransIDVec, Gap3CrossTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap3CrossTransIDBack));
    finalFlyStruct(flyStructCounter).Gap3CrossIndex = ...
        sort([strfind(CompTransIDVec, Gap3CrossTransIDFor) strfind(CompTransIDVec, Gap3CrossTransIDBack)]);
    Gap4CrossCount = length(strfind(CompTransIDVec, Gap4CrossTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap4CrossTransIDBack));
    finalFlyStruct(flyStructCounter).Gap4CrossIndex = ...
        sort([strfind(CompTransIDVec, Gap4CrossTransIDFor) strfind(CompTransIDVec, Gap4CrossTransIDBack)]);
    
    % Count number of times the TransIDs happen for gap retreats per width
    % and save the positions of the crossing events
    Gap1RetCount = length(strfind(CompTransIDVec, Gap1RetTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap1RetTransIDBack));
    finalFlyStruct(flyStructCounter).Gap1RetIndex = ...
        sort([strfind(CompTransIDVec, Gap1RetTransIDFor) strfind(CompTransIDVec, Gap1RetTransIDBack)]);
    Gap2RetCount = length(strfind(CompTransIDVec, Gap2RetTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap2RetTransIDBack));
    finalFlyStruct(flyStructCounter).Gap2RetIndex = ...
        sort([strfind(CompTransIDVec, Gap2RetTransIDFor) strfind(CompTransIDVec, Gap2RetTransIDBack)]);
    Gap3RetCount = length(strfind(CompTransIDVec, Gap3RetTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap3RetTransIDBack));
    finalFlyStruct(flyStructCounter).Gap3RetIndex = ...
        sort([strfind(CompTransIDVec, Gap3RetTransIDFor) strfind(CompTransIDVec, Gap3RetTransIDBack)]);
    Gap4RetCount = length(strfind(CompTransIDVec, Gap4RetTransIDFor)) + ...
        length(strfind(CompTransIDVec, Gap4RetTransIDBack));
    finalFlyStruct(flyStructCounter).Gap4RetIndex = ...
        sort([strfind(CompTransIDVec, Gap4RetTransIDFor) strfind(CompTransIDVec, Gap4RetTransIDBack)]);

    % Count number of times the TransIDs happen for gap circumvention per width
    % and save the positions of the crossing events
    Gap1CircCount = length(strfind(CompTransIDVec, Gap1CircTransID));
    finalFlyStruct(flyStructCounter).Gap1CircIndex = strfind(CompTransIDVec, Gap1CircTransID);
    Gap2CircCount = length(strfind(CompTransIDVec, Gap2CircTransID));
    finalFlyStruct(flyStructCounter).Gap2CircIndex = strfind(CompTransIDVec, Gap2CircTransID);
    Gap3CircCount = length(strfind(CompTransIDVec, Gap3CircTransID));
    finalFlyStruct(flyStructCounter).Gap3CircIndex = strfind(CompTransIDVec, Gap3CircTransID);
    Gap4CircCount = length(strfind(CompTransIDVec, Gap4CircTransID));
    finalFlyStruct(flyStructCounter).Gap4CircIndex = strfind(CompTransIDVec, Gap4CircTransID);
    
    % Fill the vectors with the counts
    FlyCrossCount(flyStructCounter,1) = Gap1CrossCount; 
    FlyRetCount(flyStructCounter,1) = Gap1RetCount;
    FlyCircCount(flyStructCounter,1) = Gap1CircCount;   
    
    FlyCrossCount(flyStructCounter,2) = Gap2CrossCount; 
    FlyRetCount(flyStructCounter,2) = Gap2RetCount;
    FlyCircCount(flyStructCounter,2) = Gap2CircCount;    
    
    FlyCrossCount(flyStructCounter,3) = Gap3CrossCount; 
    FlyRetCount(flyStructCounter,3) = Gap3RetCount;
    FlyCircCount(flyStructCounter,3) = Gap3CircCount;    
    
    FlyCrossCount(flyStructCounter,4) = Gap4CrossCount; 
    FlyRetCount(flyStructCounter,4) = Gap4RetCount;
    FlyCircCount(flyStructCounter,4) = Gap4CircCount;
end

% Compute the crossing rate and the binomial error per fly
FlyCrossCountRate = FlyCrossCount./(FlyRetCount+FlyCircCount+FlyCrossCount);
FlyCrossBinomErr = sqrt(FlyCrossCountRate.*(1-FlyCrossCountRate)./(FlyCrossCount+FlyRetCount+FlyCircCount));

% Compute the population mean and errors, using # of flies as n
meanCrossRate = mean(FlyCrossCountRate);
stderror = std(FlyCrossCountRate)/sqrt(length(FlyCrossCountRate));

% % Compute the population mean and errors, using # of crossings as n
% % I don't think this is the right thing to do because this likely
% % underestimates our error
% SumCrossCountRate = sum(FlyCrossCount)./(sum(FlyCrossCount)+sum(FlyRetCount)+sum(FlyCircCount));
% SumCrossBinomErr = sqrt(SumCrossCountRate.*(1-SumCrossCountRate)./(sum(FlyCrossCount)+sum(FlyRetCount)+sum(FlyCircCount)));

% Plot the population mean and error with n = # of flies
fig1 = figure(1);
errorbar(2.5:0.5:4,meanCrossRate,stderror);
xlim([2.3 4.2]);
xticks(2.5:0.5:4);
xlabel('Gap width (mm)');
ylabel('Success Rate');

% % Plot the population mean and error with n = # of crossings
% figure(2)
% errorbar(2.5:0.5:4,SumCrossCountRate,SumCrossBinomErr);
% xlim([2.3 4.2]);
% xticks(2.5:0.5:4);
% xlabel('Gap width (mm)');
% ylabel('Success Rate');

% Plot the mean and error for each fly
fig2 = figure(2);
hold on
for i = 1:length(finalFlyStruct)
errorbar(2.5:0.5:4, FlyCrossCountRate(i,:),FlyCrossBinomErr(i,:));
end
hold off
xlim([2.3 4.2]);
xticks(2.5:0.5:4);
xlabel('Gap width (mm)');
ylabel('Success Rate');

end
