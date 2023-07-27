% Choose ranges of colNum and rowNum that you want to fill
for colNum = [2, 3]
    for rowNum = [32, 40]
        % Skip all entries that don't have an associated AllCrossStatsName
        if isempty(AllCrossStatsNames{rowNum,colNum})
            pause(1); % Put in a pause to give users chance to abort
            continue
        end
        % Ask user to choose the combined WS for the appropriate genotype
        [file,path] = ...
            uigetfile('D:\Gap_Crossing\Data',...
                      ['Select ', AllCrossStatsNames{rowNum,colNum}, ' Combined WS']);
        % If the user skips, move on to the next iteration of the loop
        if file == 0
            pause(1); % Put in a pause to give users chance to abort
            continue
        end
        % Load in the WS, strip the extra layer, then save crossStats and
        % genotype to AllCrossStats and AllCrossStatsGenotypes
        tempWS = load(fullfile(path,file));
        tempWS = tempWS.WS;
        AllCrossStats{rowNum,colNum} = tempWS.crossStats;
        AllCrossStatsGenotypes{rowNum,colNum} = tempWS.genotype;
        % Output a message with which file was loaded so users can double
        % check that they loaded in what they thought they did
        disp(['Loaded file: ', file]);
    end
end