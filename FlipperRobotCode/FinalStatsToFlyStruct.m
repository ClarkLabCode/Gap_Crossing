% Converts the finalStats struct into finalFlyStruct which is a fly-centric
% structure

function finalFlyStruct = FinalStatsToFlyStruct(finalStats, NumCorridors)

% Sort finalStats by CorridorID (needs to be done in table format first)
finalStatsTab = struct2table(finalStats);
finalStatsTab = sortrows(finalStatsTab,'CorridorID');
finalStatsNew = table2struct(finalStatsTab);

% Delete rows where flies are not in any of the corridors
CorrVec = [finalStatsNew(:).CorridorID];
CorrVecLog = (CorrVec == 0);
finalStatsNew(CorrVecLog) = [];

% Remove any rows which contain more than 1 fly in a particular corridor
% during any one frame of the video
IDVec = [finalStatsNew(:).CorridorID; finalStatsNew(:).FlipNumber; finalStatsNew(:).FrameInFlip]';
[~, ~, ic] = unique(IDVec, 'rows', 'stable');
[count, ~] = hist(ic, unique(ic));
CountLog = ismember(ic,find(count==1));
finalStatsNew = finalStatsNew(CountLog);

% Get the necessary info to know which rows correspond to which frame
FinalCorrVec = [finalStatsNew(:).CorridorID];
[numRowsPerCorr, ~] = hist(FinalCorrVec, 1:NumCorridors);

% Define finalFlyStruct
finalFlyStruct = struct();

% Used for indexing through finalStats
LowerLim = 1;
UpperLim = numRowsPerCorr(1);

% Go through finalStats and transfer the data into finalFlyStruct
for CorrCounter = 1:NumCorridors
    ii = 1;
    for Row = LowerLim:UpperLim
        finalFlyStruct(CorrCounter).Area(ii) = finalStatsNew(Row).Area;
        finalFlyStruct(CorrCounter).CentroidX(ii) = finalStatsNew(Row).Centroid(1);
        finalFlyStruct(CorrCounter).CentroidY(ii) = finalStatsNew(Row).Centroid(2);
        finalFlyStruct(CorrCounter).AbsoluteTime(ii) = finalStatsNew(Row).AbsoluteTime;
        finalFlyStruct(CorrCounter).BoundingBoxTLX(ii) = finalStatsNew(Row).BoundingBox(1);
        finalFlyStruct(CorrCounter).BoundingBoxTLY(ii) = finalStatsNew(Row).BoundingBox(2);
        finalFlyStruct(CorrCounter).BoundingBoxWidth(ii) = finalStatsNew(Row).BoundingBox(3);
        finalFlyStruct(CorrCounter).BoundingBoxHeight(ii) = finalStatsNew(Row).BoundingBox(4);
        finalFlyStruct(CorrCounter).MajorAxisLength(ii) = finalStatsNew(Row).MajorAxisLength;
        finalFlyStruct(CorrCounter).MinorAxisLength(ii) = finalStatsNew(Row).MinorAxisLength;
        finalFlyStruct(CorrCounter).Orientation(ii) = finalStatsNew(Row).Orientation;
        finalFlyStruct(CorrCounter).FlipNumber(ii) = finalStatsNew(Row).FlipNumber;
        finalFlyStruct(CorrCounter).FrameInFlip(ii) = finalStatsNew(Row).FrameInFlip;
        finalFlyStruct(CorrCounter).CompID(ii) = finalStatsNew(Row).CompID;
        ii = ii + 1;
    end
    % LowerLim and UpperLim are how we find which index we need from finalStats
    LowerLim = LowerLim+numRowsPerCorr(CorrCounter);
    if CorrCounter ~= NumCorridors
        UpperLim = UpperLim+numRowsPerCorr(CorrCounter+1);
    end
    % finalFlyStruct(CorrCounter).CassetteID = finalStatsNew(Row).CassetteID;
    finalFlyStruct(CorrCounter).CorridorID = finalStatsNew(Row).CorridorID;
end

end
