function GapCrossingAnalysis()

% dlgList contains a list of all the input parameters to be asked of user
dlgList = {'Genotype','Date','Time','Temperature'}

% genotype = 'IsoD1';
% dateAcq = '12-29-20';
% timeAcq = '16-35-12';
% eclosionDate = '12-25-20';
% timeZone = 15;
% notes = ' ';
% flipRate = 30 (in seconds);
% expLength = 60 (in mins);
% directoryName = 'test';
% sizeThreshCutOff = 100;
% indPosFrameBuffer = 5;

expParams = inputdlg(dlgList);

% Need to convert numeric entries using str2num()

% Runs all gap crossing scripts in appropriate order