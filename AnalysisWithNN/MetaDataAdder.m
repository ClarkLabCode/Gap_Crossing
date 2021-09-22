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

function finalFlyStruct = MetaDataAdder(finalFlyStruct, genotype, dateAcq,...
    timeAcq, eclosionDate, timeZone, notes, flipRate, expLength, temp, CassetteID, ...
    directoryName, sizeThreshCutOff, indPosFrameBuffer, ExpRepNum, inputFileName)

for flyStructCounter = 1:length(finalFlyStruct)
    finalFlyStruct(flyStructCounter).Genotype = genotype;
    finalFlyStruct(flyStructCounter).DateAcq = dateAcq;
    finalFlyStruct(flyStructCounter).TimeAcq = timeAcq;
    finalFlyStruct(flyStructCounter).EclosionDate = eclosionDate;
    finalFlyStruct(flyStructCounter).TimeZone = timeZone;
    finalFlyStruct(flyStructCounter).Notes = notes;
    finalFlyStruct(flyStructCounter).FlipRate = flipRate;
    finalFlyStruct(flyStructCounter).ExpLength = expLength;
    finalFlyStruct(flyStructCounter).Temperature = temp;
    finalFlyStruct(flyStructCounter).CassetteID = CassetteID;
    finalFlyStruct(flyStructCounter).DirectoryName = directoryName;
    finalFlyStruct(flyStructCounter).SizeThreshCutOff = sizeThreshCutOff;
    finalFlyStruct(flyStructCounter).IndPosFrameBuffer = indPosFrameBuffer;
    finalFlyStruct(flyStructCounter).ExpNum = ExpRepNum;
    finalFlyStruct(flyStructCounter).VidFile = inputFileName;
end
