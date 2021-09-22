% Removes flies from finalFlyStruct that do not move for at least 20% of
% the video

function finalFlyStruct = FlyActivityFilter(finalFlyStruct)

maxTimePts = 0;

% Figure out how many frames the most active fly was active for
for flyStructCounter = 1:length(finalFlyStruct)
    maxTimePts = max(size(finalFlyStruct(flyStructCounter).Area,2),maxTimePts);
end

timePtsFilter = zeros(1,length(finalFlyStruct));

% Determine which flies were active for at least 20% of the time
for flyStructCounter = 1:length(finalFlyStruct)
    timePtsFilter(flyStructCounter) = size(finalFlyStruct(flyStructCounter).Area,2) > maxTimePts/5;
end

% Convert this to logical
timePtsFilter = logical(timePtsFilter);

% Only keep the flies that were active for at least 20%
finalFlyStruct = finalFlyStruct(timePtsFilter);

end