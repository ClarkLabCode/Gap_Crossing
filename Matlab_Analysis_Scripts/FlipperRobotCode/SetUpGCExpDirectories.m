% SETUPGCEXPDIRECTORIES Sets up all necessary subdirectories for saving data
%
%  Makes the necessary folders for progress to be saved at each checkpoint
%  in FullGapCrossingAnalysis. These folders are:
%    0_Raw_Videos
%    1_Params_Entered
%    2_Video_Processed
%    3_Corridor_Skeleton_Made
%    4_Neural_Net_Analyzed
%    5_Fully_Analyzed
%    Background_Frames
%    Fly_Structure
%    Plots


function WS_names = SetUpGCExpDirectories(WS_all,LocalGCDirectoryPath)

% Establish number of GC experiments
NumGCExp = length(WS_all);

% Initialize WS_names to hold all the WS file names that are being made
WS_names = cell(NumGCExp,1);

% Loop through all experiments loaded in
for vidCounter = 1:NumGCExp

    % Navigate to the Data folder
    cd([LocalGCDirectoryPath,'Data'])

    % Load in WS for the appropriate experiment
    WS = WS_all{vidCounter};

    % Make folder with Directory Name
    mkdir(WS.directoryName);

    % Navigate into new folder
    cd(WS.directoryName);

    % Make folder with experiment number
    mkdir(['Experiment_', num2str(WS.ExpRepNum)]);

    % Navigate into new folder
    cd(['Experiment_', num2str(WS.ExpRepNum)]);

    % Make the appropriate subfolders
    mkdir '0_Raw_Videos';
    mkdir '1_Params_Entered';
    mkdir '2_Video_Processed';
    mkdir '3_Corridor_Skeleton_Made';
    mkdir '4_Neural_Net_Analyzed';
    mkdir '5_Fully_Analyzed';
    mkdir 'Background_Frames';
    mkdir 'Fly_Structure';
    mkdir 'Plots';

    % Navigate into Raw_Videos and copy in the video file
    cd '0_Raw_Videos';
    copyfile([LocalGCDirectoryPath,'Data\All_Raw_Videos\', WS.inputFileName]);

    % Navigate into Params_Entered and save WS
    cd ../1_Params_Entered
    save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS')

    % Now fill WS_names appropriately
    WS_names{vidCounter} = ...
        ['Data\',WS.directoryName, ...
         '\Experiment_', num2str(WS.ExpRepNum),'\1_Params_Entered\', ...
         WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'];

    % Now clear WS
    clear WS
    
    % Navigate back into the base code directory
    cd([LocalGCDirectoryPath,'Matlab_Analysis_Scripts\FlipperRobotCode']);

end

end