% FLYACTIVITYFILTER Removes inactive fly statistics from finalFlyStruct
%
%  Removes flies from finalFlyStruct that do not move for at least 20% of
%  the video. This is to prevent inactive flies with little data from
%  potentially disproportionately affecting population crossing statistics.

function WS = FlyActivityFilter(WS)

% function finalFlyStruct = FlyActivityFilter(finalFlyStruct)

% Port in the relevant fields from WS
finalFlyStruct = WS.finalFlyStruct;

maxTimePts = 0;

% Figure out how many frames the most active fly was active for
for flyStructCounter = 1:length(finalFlyStruct)
    maxTimePts = max(size(finalFlyStruct(flyStructCounter).Area,2),maxTimePts);
end

% Make sure that maxTimePts is at least 1/2 of the length of the video
maxTimePts = max(maxTimePts, round(length(WS.frameMarkerVec)/2));

% Create an array that will hold whether or not a fly is active enough
timePtsFilter = zeros(1,length(finalFlyStruct));

% Determine which flies were active for at least 20% of the time
for flyStructCounter = 1:length(finalFlyStruct)
    timePtsFilter(flyStructCounter) = size(finalFlyStruct(flyStructCounter).Area,2) > maxTimePts/5;
end

% Convert this to logical
timePtsFilter = logical(timePtsFilter);

% Only keep the flies that were active for at least 20%
finalFlyStruct = finalFlyStruct(timePtsFilter);

% Update the fields in WS
WS.finalFlyStruct = finalFlyStruct;

end
