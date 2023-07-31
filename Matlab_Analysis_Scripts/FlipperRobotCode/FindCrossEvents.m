% FINDCROSSEVENTS Identifies all crossing events and adds fields to FlipBinnedFlyStruct
%
%  Finds all crossings (glass and proper are combined at this stage of
%  analysis until the NN is used to classify the two), retreats, and 
%  circumventions per fly per gap width and adds them to FlipBinnedFlyStruct 
%  in the form of new fields. This is done by looking for the appropriate
%  compartment transitions for each fly. For a reminder of how a compartment 
%  is defined, please refer to the diagram in the README.

function WS = FindCrossEvents(WS)

% function FlipBinnedFlyStruct = FindCrossEvents(FlipBinnedFlyStruct, NumGaps)

% Port in the relevant fields from WS
FlipBinnedFlyStruct     = WS.FlipBinnedFlyStruct;
NumGaps                 = WS.NumGaps;

% Compute number of unique compartments
NumComps = 3*NumGaps+1;

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

% % % % % Create vectors to hold cross, ret, circ info for each fly per gap width
% % % % FlyCrossCount = zeros(length(FlipBinnedFlyStruct.ExpNum),NumGaps);
% % % % FlyRetCount = zeros(length(FlipBinnedFlyStruct.ExpNum),NumGaps);
% % % % FlyCircCount = zeros(length(FlipBinnedFlyStruct.ExpNum),NumGaps);

for flyStructCounter = 1:length(FlipBinnedFlyStruct.ExpNum)
    for OddFlipCounter = 1:length(FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips)
        % Find the frames in which a fly changes compartments
        CompTransVec = [1 diff([FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).CompID])];
        % Turn these into logical arrays
        CompTransLog = CompTransVec ~= 0;
        % Declare new field for UniqCompID and give placeholder value
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UniqCompID = 0;
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UniqCompIDIndex = 1; %Placeholder 1 since index
        % Only grab unique CompIDs and save those and their index
        % Check to make sure nothing is empty first
        if sum(CompTransLog) ~= 0
            CompTransIDVec = FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).CompID(CompTransLog);
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UniqCompID = CompTransIDVec;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UniqCompIDIndex = find(CompTransLog);
        end
        % Declare new fields for holding # of crosses, rets, circs
        % Give placeholder values of 0
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrosses = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreats = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventions = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrosses = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreats = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventions = zeros(1,NumGaps);
        % Declare new fields for holding index of locations of gap events
        % Give placeholder values of 0
        for GapNumber = 1:NumGaps
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrossesIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreatsIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventionsIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrossesIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreatsIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventionsIndex(GapNumber).GapID = 0;
        end
        % Count number of crosses, retreats, and circumventions in that flip
        % and mark the index locations (if they occur)
        for GapNumber = 1:NumGaps
            if strcmp(FlipBinnedFlyStruct.ExpNum(flyStructCounter).IdData.LocOfSmallestGapInOddFlip, 'Bottom')
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,1)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,1)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,2)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,1)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,1)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,2)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,1)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,2)');
                end
            else
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,2)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,2)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,1)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,2)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,2)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,1)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).UpCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,2)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.OddFlips(OddFlipCounter).DownCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,1)');
                end
            end
        end
    end
    for EvenFlipCounter = 1:length(FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips)
        % Find the frames in which a fly changes compartments
        CompTransVec = [1 diff([FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).CompID])];
        % Turn these into logical arrays
        CompTransLog = CompTransVec ~= 0;
        % Declare new field for UniqCompID and give placeholder value
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UniqCompID = 0;
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UniqCompIDIndex = 1; %Placeholder 1 since index
        % Only grab unique CompIDs and save those and their index
        % Check to make sure nothing is empty first
        if sum(CompTransLog) ~= 0
            CompTransIDVec = FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).CompID(CompTransLog);
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UniqCompID = CompTransIDVec;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UniqCompIDIndex = find(CompTransLog);
        end
        % Declare new fields for holding # of crosses, rets, circs 
        % Give placeholder values of 0
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrosses = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreats = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventions = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrosses = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreats = zeros(1,NumGaps);
        FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventions = zeros(1,NumGaps);
        % Declare new fields for holding index of locations of gap events
        % Give placeholder values of 0
        for GapNumber = 1:NumGaps
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrossesIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreatsIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventionsIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrossesIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreatsIndex(GapNumber).GapID = 0;
            FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventionsIndex(GapNumber).GapID = 0;
        end
        % Count number of crosses, retreats, and circumventions in that flip
        % and mark the index locations (if they occur)
        for GapNumber = 1:NumGaps
            if strcmp(FlipBinnedFlyStruct.ExpNum(flyStructCounter).IdData.LocOfSmallestGapInOddFlip, 'Top')
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,1)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,1)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,2)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,1)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,1)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,2)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,1)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,2)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,1)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,2)');
                end
            else
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,2)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrosses(GapNumber) = ...
                    length(strfind(CompTransIDVec, CrossID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,2)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrosses(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCrossesIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CrossID(:,GapNumber,1)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,2)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreats(GapNumber) = ...
                    length(strfind(CompTransIDVec, RetID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,2)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreats(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownRetreatsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, RetID(:,GapNumber,1)');
                end
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,2)'));
                FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventions(GapNumber) = ...
                    length(strfind(CompTransIDVec, CircID(:,GapNumber,1)'));
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).UpCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,2)');
                end
                if FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventions(GapNumber) ~= 0
                    FlipBinnedFlyStruct.ExpNum(flyStructCounter).BehavData.EvenFlips(EvenFlipCounter).DownCircumventionsIndex(GapNumber).GapID = ...
                        strfind(CompTransIDVec, CircID(:,GapNumber,1)');
                end
            end
        end
    end
                   
end

% Update the fields in WS
WS.FlipBinnedFlyStruct = FlipBinnedFlyStruct;

end
