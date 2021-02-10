% Takes the finalFlyStruct and rearranges the data such that it separates
% out behavior data between flips to make it easier to manage

function NewFlyStruct = FlipBinner(finalFlyStruct)

for flyStructCounter = 1:length(finalFlyStruct)
    
    % NewFlyStruct.Metadata
    NewFlyStruct(flyStructCounter).MetaData.DirectoryName = ...
        finalFlyStruct(flyStructCounter).DirectoryName;
    NewFlyStruct(flyStructCounter).MetaData.Genotype = ...
        finalFlyStruct(flyStructCounter).Genotype;
    NewFlyStruct(flyStructCounter).MetaData.DateAcq = ...
        finalFlyStruct(flyStructCounter).DateAcq;
    NewFlyStruct(flyStructCounter).MetaData.TimeAcq = ...
        finalFlyStruct(flyStructCounter).TimeAcq;
    NewFlyStruct(flyStructCounter).MetaData.EclosionDate = ...
        finalFlyStruct(flyStructCounter).EclosionDate;
    NewFlyStruct(flyStructCounter).MetaData.TimeZone = ...
        finalFlyStruct(flyStructCounter).TimeZone;
    NewFlyStruct(flyStructCounter).MetaData.FlipRate = ...
        finalFlyStruct(flyStructCounter).FlipRate;
    NewFlyStruct(flyStructCounter).MetaData.ExpLength = ...
        finalFlyStruct(flyStructCounter).ExpLength;
    NewFlyStruct(flyStructCounter).MetaData.SizeThreshCutOff = ...
        finalFlyStruct(flyStructCounter).SizeThreshCutOff;
    NewFlyStruct(flyStructCounter).MetaData.IndPosFrameBuffer = ...
        finalFlyStruct(flyStructCounter).IndPosFrameBuffer;
    NewFlyStruct(flyStructCounter).MetaData.Notes = ...
        finalFlyStruct(flyStructCounter).Notes;
    
    % NewFlyStruct.IdData
    NewFlyStruct(flyStructCounter).IdData.CassetteID = ...
        finalFlyStruct(flyStructCounter).CassetteID;
    NewFlyStruct(flyStructCounter).IdData.CorridorID = ...
        finalFlyStruct(flyStructCounter).CorridorID;
    
    % NewFlyStruct.BehavData
    for FlipCounter = 1:max(finalFlyStruct(flyStructCounter).FlipNumber)
        FlipLog = (finalFlyStruct(flyStructCounter).FlipNumber == FlipCounter);
        % This is where I left off
    end
    
end

end