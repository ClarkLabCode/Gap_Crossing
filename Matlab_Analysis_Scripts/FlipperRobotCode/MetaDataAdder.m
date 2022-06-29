% Adds meta data info to the struct

% Example inputs:

% genotype = 'IsoD1';
% dateAcq = '12-29-20';
% timeAcq = '16-35-12';
% eclosionDate = '12-25-20';
% timeZone = 15;
% notes = ' ';
% flipRate = 15 (in seconds);
% expLength = 60 (in mins);
% temp = 30 (in C);
% CassetteID = 1;
% directoryName = 'IsoD1_BC';
% sizeThreshCutOff = 100;
% indPosFrameBuffer = 5;
% ExpRepNum = 1;
% inputFileName = '2020-12-29 16-35-22_IsoD1_30C_MF.mov';

function WS = MetaDataAdder(WS)

% function finalFlyStruct = MetaDataAdder(finalFlyStruct, genotype, dateAcq,...
%     timeAcq, eclosionDate, timeZone, notes, flipRate, expLength, temp, CassetteID, ...
%     directoryName, sizeThreshCutOff, indPosFrameBuffer, ExpRepNum, inputFileName, ...
%     LocOfSmallGapInOddFlip)

% Port in the relevant fields from WS
finalFlyStruct = WS.finalFlyStruct;
genotype = WS.genotype;
dateAcq = WS.dateAcq;
timeAcq = WS.timeAcq;
eclosionDate = WS.eclosionDate;
timeZone = WS.timeZone;
notes = WS.notes;
flipRate = WS.flipRate;
expLength = WS.expLength;
temperature = WS.temperature;
directoryName = WS.directoryName;
sizeThreshCutOff = WS.sizeThreshCutOff;
indPosFrameBuffer = WS.indPosFrameBuffer;
ExpRepNum = WS.ExpRepNum;
inputFileName = WS.inputFileName;
LocOfSmallestGapInOddFlip = WS.LocOfSmallestGapInOddFlip;

% Now fill them into the fly struct
for flyStructCounter = 1:length(finalFlyStruct)
    finalFlyStruct(flyStructCounter).Genotype = genotype;
    finalFlyStruct(flyStructCounter).DateAcq = dateAcq;
    finalFlyStruct(flyStructCounter).TimeAcq = timeAcq;
    finalFlyStruct(flyStructCounter).EclosionDate = eclosionDate;
    finalFlyStruct(flyStructCounter).TimeZone = timeZone;
    finalFlyStruct(flyStructCounter).Notes = notes;
    finalFlyStruct(flyStructCounter).FlipRate = flipRate;
    finalFlyStruct(flyStructCounter).ExpLength = expLength;
    finalFlyStruct(flyStructCounter).Temperature = temperature;
    finalFlyStruct(flyStructCounter).DirectoryName = directoryName;
    finalFlyStruct(flyStructCounter).SizeThreshCutOff = sizeThreshCutOff;
    finalFlyStruct(flyStructCounter).IndPosFrameBuffer = indPosFrameBuffer;
    finalFlyStruct(flyStructCounter).ExpNum = ExpRepNum;
    finalFlyStruct(flyStructCounter).VidFile = inputFileName;
    finalFlyStruct(flyStructCounter).LocOfSmallestGapInOddFlip = LocOfSmallestGapInOddFlip;
end

% Update the fields in WS
WS.finalFlyStruct = finalFlyStruct;

end