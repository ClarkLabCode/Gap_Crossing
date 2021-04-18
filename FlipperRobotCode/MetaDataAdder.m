% Adds meta data info to the struct

% Example inputs:

% genotype = 'IsoD1';
% dateAcq = '12-29-20';
% timeAcq = '16-35-12';
% eclosionDate = '12-25-20';
% timeZone = 15;
% notes = ' ';
% flipRate = 30 (in seconds);
% expLength = 60 (in mins);
% temp = 30 (in C);
% directoryName = 'test';
% sizeThreshCutOff = 100;
% indPosFrameBuffer = 5;

function finalFlyStruct = MetaDataAdder(finalFlyStruct, genotype, dateAcq,...
    timeAcq, eclosionDate, timeZone, notes, flipRate, expLength, temp, ...
    directoryName, sizeThreshCutOff, indPosFrameBuffer)

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
    finalFlyStruct(flyStructCounter).DirectoryName = directoryName;
    finalFlyStruct(flyStructCounter).SizeThreshCutOff = sizeThreshCutOff;
    finalFlyStruct(flyStructCounter).IndPosFrameBuffer = indPosFrameBuffer;
end
