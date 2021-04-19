% Runs all gap crossing scripts in appropriate order and saves all data
% into the subfolder structure: Data\[Directory Name]\...

% VidPrevProcessed = Whether the raw video has been processed before
% Since this is by far the most time-consuming part of the analysis
% pipeline, this option is here to allow you to skip it if you ran it
% separately

% VidPrevProcessed = 1 skips FindFramesForFlip.m & RawFlyInfoExtracter.m
% and then requests the user provides the location of the analyzed videos

% Use the below command in the command window to navigate to this script
% cd 'C:\Users\clarklab\Joe\Gap_Crossing\Matlab_Analysis_Scripts\FlipperRobotCode'

function GapCrossingAnalysis(VidPrevProcessed)

disp('Please make sure the workspace is empty before proceeding.');
disp('If it is not empty, cancel the next two pop-up boxes and enter "clear" in the command window.');

% Add path to directory with all analysis scripts and navigate to it
addpath 'C:\Users\clarklab\Joe\Gap_Crossing\Matlab_Analysis_Scripts\FlipperRobotCode'
cd 'C:\Users\clarklab\Joe\Gap_Crossing\Matlab_Analysis_Scripts\FlipperRobotCode'

% Create function handles for all necessary functions so that they are
% accessible outside of the analysis script directory
FindFramesForFlipHan = @FindFramesForFlip;
RawFlyInfoExtracterHan = @RawFlyInfoExtracter;
CorridorIdentifierHan = @CorridorIdentifier;
CompTracerHan = @CompTracer;
CompIdentifierHan = @CompIdentifier;
FinalStatsToFlyStructHan = @FinalStatsToFlyStruct;
FlyActivityFilterHan = @FlyActivityFilter;
FindCrossEventsAndStatsHan = @FindCrossEventsAndStats;
MetaDataAdderHan = @MetaDataAdder;
% FlipBinnerHan = @FlipBinner;

% Navigate to directory with the raw videos to be processed
cd ..\..\Data\All_Raw_Videos

[file, path] = uigetfile('*.*');
absFilePath = fullfile(path,file);
inputFileName = absFilePath;

% dlgList contains a list of all the input parameters to be asked of user
dlgList = {'Genotype','Directory Name','Input File Name',...
    'Flip Rate (sec)','Experiment Length (min)','Temperature (C)',...
    'Date Acquired','Time Acquired','Eclosion Date','Time Zone','Notes',...
    'Cassette ID','Number of Corridors','Number of Gaps per Corridor',...
    'sizeThreshCutOff','indPosFrameBuffer','Number of Pixels to Erode'};
dlgtitle = 'Input';

% Set up dimensions of inputdlg box
dimsH = [1,1,4,...
    1,1,1,...
    1,1,1,1,5,...
    1,1,1,...
    1,1,1];
dimsW = ones(1,length(dimsH))*50;
dims = [dimsH;dimsW]';

% Define default inputs for inputdlg box
definput = {'','',inputFileName,...
    '30','60','30',...
    '','','','','',...
    '1','7','4',...
    '100','5','1'};

% Pop up a dialog box to fill out all info into the fields
expParams = inputdlg(dlgList,dlgtitle,dims,definput,'on');

% Convert inputs into variables
genotype            = expParams{1};
directoryName       = expParams{2};
inputFileName       = expParams{3};
flipRate            = str2double(expParams{4});
expLength           = str2double(expParams{5});
temp                = str2double(expParams{6});
dateAcq             = expParams{7};
timeAcq             = expParams{8};
eclosionDate        = expParams{9};
timeZone            = expParams{10};
notes               = expParams{11};
CassetteID          = str2double(expParams{12});
NumCorridors        = str2double(expParams{13});
NumGaps             = str2double(expParams{14});
sizeThreshCutOff    = str2double(expParams{15});
indPosFrameBuffer   = str2double(expParams{16});
erodePix            = str2double(expParams{17});

% Checks to make sure that no critical fields are empty
if (isempty(expParams{1}) || isempty(expParams{2}) || isempty(expParams{3}))
    error('Genotype, Directory Name, and/or File Name field(s) empty.');
end

% Navigate up one level from directory with the raw videos
cd ..

% Make folder with Directory Name
mkdir(directoryName);

% Navigate into new folder
cd(directoryName);

% Make the appropriate subfolders
mkdir 'Background_Frames';
mkdir 'Fly_Structure';
mkdir 'Processed_Video_Structure';
mkdir 'Raw_Videos';

% Navigate into Raw_Videos and copy in the video file
cd 'Raw_Videos';
copyfile(inputFileName);

% If video was previously processed, load in the workspace
if VidPrevProcessed == 1
    cd ../Processed_Video_Structure
    [procFile,procPath] = uigetfile('*.*','Select Processed Workspace File');
    absProcFilePath = fullfile(procPath,procFile);
    load(absProcFilePath);
% Else, process the video and save the processed workspace
else
% Navigate into Background_Frames and runs video processing scripts
% (Video processing = FindFramesForFlip.m & RawFlyInfoExtracter.m)
% After these run, the raw video is no longer necessary and all info is
% stored in finalStats structure which later gets converted to finalFlyStruct
cd ../Background_Frames
[indPos, ~, frameMarkerVec, changeVecDil, dotProdVec] = ...
    FindFramesForFlipHan(inputFileName,flipRate,0);
[finalStats, AreaVec, AreaVecLog] = ...
    RawFlyInfoExtracterHan(inputFileName,directoryName,...
    sizeThreshCutOff,indPosFrameBuffer,erodePix,indPos);

% Navigate to Processed_Video_Structure and save finalStats
cd ../Processed_Video_Structure
save([directoryName,'_finalStats.mat'],'finalStats');
save([directoryName,'_workspace.mat']);

end

% Now run the rest of the analysis
[finalStats, CorrMask1, CorrMask2] = ...
    CorridorIdentifierHan(inputFileName, finalStats, NumCorridors, indPos);
[CompMask1, CompMask2, NumComps] = CompTracerHan(inputFileName, NumCorridors, NumGaps, indPos);
finalStats = CompIdentifierHan(finalStats, CompMask1, CompMask2, indPos, NumComps);
finalFlyStruct = FinalStatsToFlyStructHan(finalStats, NumCorridors);
finalFlyStruct = FlyActivityFilterHan(finalFlyStruct);
[finalFlyStruct, meanCrossRate, stderror, FlyCrossCountRate, FlyCrossBinomErr] = ...
    FindCrossEventsAndStatsHan(finalFlyStruct, NumGaps);
finalFlyStruct = MetaDataAdderHan(finalFlyStruct, genotype, dateAcq,...
    timeAcq, eclosionDate, timeZone, notes, flipRate, expLength, temp, CassetteID, ...
    directoryName, sizeThreshCutOff, indPosFrameBuffer);
% FlipBinnerHan()

% Navigate to Fly_Structure and save structure there
cd ../Fly_Structure
save([directoryName,'_finalFlyStruct.mat'],'finalFlyStruct');
save([directoryName,'_analyzed_workspace.mat']);

end
