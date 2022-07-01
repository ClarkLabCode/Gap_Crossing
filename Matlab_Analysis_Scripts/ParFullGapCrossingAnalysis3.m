function ParFullGapCrossingAnalysis3(varargin)

% Parse the input to see if user specified how many experiments are to be
% analyzed or not, and throw an error if multiple inputs were given or if
% the input given was not a number
if nargin > 1
    error('Unexpected number of inputs.')
elseif nargin == 1
    if ~isa(varargin{1},'double')
        error(['Expected function argument to be empty ', ...
               'or equal to the number of experiments to be analyzed.'])
    end
end

% Create function handles for all the functions being used so that they're
% easy to identify and find for debugging reasons and so they're globally
% available within this function call
% To search for all nested functions in this code, use Ctrl+F and search "Func"
GCADirectoryCheckerFunc             = @GCADirectoryChecker;
LoadGCExpWSNamesFunc                = @LoadGCExpWSNames;
RequestUserGCInputsFunc             = @RequestUserGCInputs;
SetUpGCExpDirectoriesFunc           = @SetUpGCExpDirectories;
SaveCheckpointPromptFunc            = @SaveCheckpointPrompt;
AnalysisStepCheckFunc               = @AnalysisStepCheck;
FindFramesForFlipFunc               = @FindFramesForFlip2;
RawFlyInfoExtracterFunc             = @RawFlyInfoExtracter;
CorridorSkeletonFinderFunc          = @CorridorSkeletonFinder;
SkeletonToCorrMasksFunc             = @SkeletonToCorrMasks;
SkeletonToGapCentersFunc            = @SkeletonToGapCenters;
SkeletonToCompMasksFunc             = @SkeletonToCompMasks;
FlyCorrIdentifierFunc               = @FlyCorrIdentifier;
FlyCompIdentifierFunc               = @FlyCompIdentifier;
FinalStatsToFlyStructFunc           = @FinalStatsToFlyStruct;
FlyActivityFilterFunc               = @FlyActivityFilter;
MetaDataAdderFunc                   = @MetaDataAdder;
FlipBinnerFunc                      = @FlipBinner2;
ZapStationaryFlipsFunc              = @ZapStationaryFlips2;
FindCrossEventsFunc                 = @FindCrossEvents;
CrossingEventFrameFinderFunc        = @CrossingEventFrameFinder;
PredictOnOffGlassFunc               = @PredictOnOffGlass;
FinalFlipRemoverFunc                = @FinalFlipRemover;
ClassifyOnOffFunc                   = @ClassifyOnOff;

% Check to make sure that we are in the correct directory to run this
% function, and if we aren't, have the user navigate to the right directory
% This also adds the path to this directory upon success of the function
LocalGCDirectoryPath = GCADirectoryCheckerFunc;

% Warn user that this function clears workspace and figures and confirm
% that this is okay before proceeding
closeCheck = ...
    questdlg('This function will close all figures. Proceed?',...
             'Close figure check', 'Yes', 'No', 'Yes');
if ~strcmp(closeCheck,'Yes')
    error('User terminated function call.')
end

% Close all figures
close all;

% Ask user how many experiments will be analyzed if it wasn't given as an
% input to the function in the first place
if nargin < 1
    NumGCExp = input('How many experiments do you want to analyze?\n');
else
    NumGCExp = varargin{1};
end

% Check to make sure that not too many experiments are being loaded that
% would overload the RAM and also that the user intended to input such a
% large number of experiments
if NumGCExp > 66
    error('Too many experiments attempting to be loaded. Can load up to 66.')
elseif NumGCExp > 20
    VerifyNumGCExp = ...
        questdlg('Are you sure you want to analyze this many experiments?',...
             'Analyzing many experiments', 'Yes', 'No', 'Yes');
    if ~strcmp(VerifyNumGCExp,'Yes')
        error(['User canceled out of function because ', ...
               'of a misinput in the number of experiments to load.'])
    end
end

% Create the options needed for save and stop features later
contOpt = 'Yes, proceed with analysis';
stopOpt = 'No, exit analysis';
SaveCheckpointResponse = contOpt;

% Create the options for analysis step detection features later
ProcessVideo = 'Process Video';
MakeSkeleton = 'Make Skeleton';
RunNeuralNet = 'Run Neural Net';
RunFinalStep = 'Run Final Steps';
AnalysisDone = 'None';

% Ask user if the experiments have previously been analyzed
prevAnalysis = ...
    questdlg('Do you want to search for and load prior analysis progress?', ...
              'Load analysis', 'Yes', 'No', 'Yes');
          
% Ask the user if they would like not to be asked to stop at checkpoints
% Note that this does not turn off intermediate saving, just the dialogue
% asking the user whether or not they would like to stop after a checkpoint
availCheckpoints = ...
    questdlg('Do you want to be offered checkpoints throughout the analysis?', ...
             'Checkpoint availability', 'Yes', 'No', 'Yes');

% If they have been previously analyzed, let user load in the WS structs
if strcmp(prevAnalysis,'Yes')
    WS_names = LoadGCExpWSNamesFunc(NumGCExp,LocalGCDirectoryPath);
% If the experiments are not previously analyzed, start requesting the
% user to fill in the necessary info
else
    % Ask user if running standard analysis 
    % (i.e., default values for: NN decision threshold, fly size threshold,
    % cassette flip frame buffer, blob detection pixel erosion size)
    StandardAnalysis = ...
        questdlg('Run analysis with default analysis values?',...
                 'Default values', 'Yes', 'No', 'Yes');

    % Convert StandardAnalysis to logical
    if strcmp(StandardAnalysis,'Yes')
        StandardAnalysis = true;
    else
        StandardAnalysis = false;
    end
    
    % Call function to request user inputs for GC analysis
    WS_all = RequestUserGCInputsFunc(NumGCExp, StandardAnalysis,LocalGCDirectoryPath);
    
    % Call function to set up the directories for all the experiments
    WS_names = SetUpGCExpDirectoriesFunc(WS_all,LocalGCDirectoryPath);
    
    % Remove WS_all since it's not needed anymore and occupies memory
    clear WS_all
    
    % Give user choice to stop here since everything has now been saved,
    % but only do this if the user indicated that they wanted checkpoints
    if ~strcmp(availCheckpoints,'No')
        SaveCheckpointResponse = SaveCheckpointPromptFunc(contOpt, stopOpt);
    end

end

% Create a loop that performs the subsequent analysis step until the user 
% requests to stop at a save checkpoint or until the files are fully analyzed
while strcmp(SaveCheckpointResponse, contOpt)
    
tic
    
    % Call a function to check what the current and next step of analysis is
    [AnalysisStep, Next_WS_names, SkeletonMade] = AnalysisStepCheckFunc(WS_names);

    % Now execute the next analysis step based on AnalysisStep
    switch AnalysisStep
        case ProcessVideo
            % Loop through all experiments loaded in
            for vidCounter = 1:NumGCExp
                % Load in WS for each experiment
                WS = load(WS_names{vidCounter});
                % To get the correct WS, we must go in one layer
                WS = WS.WS;

                % Move to the appropriate folder for this experiment
                cd([LocalGCDirectoryPath,'Data\',...
                    WS.directoryName,'\Experiment_', num2str(WS.ExpRepNum)])

                % Navigate to Background_Frames so the bg frames are saved there
                cd Background_Frames

                % Run FindFramesForFlip and RawFlyInfoExtracter
                WS = FindFramesForFlipFunc(WS);
                WS = RawFlyInfoExtracterFunc(WS);

                % If skeleton was made before video was processed, we want to
                % save WS now in Skeleton_Made because it has now had both the
                % video processed and the skeleton made, but we then also want
                % to save WS in Video_Processed without the skeleton info so
                % that we can give ourselves the option to remake skeleton
                % later if we choose to do so without rerunning video processing
                if SkeletonMade
                    % Navigate to Skeleton_Made and save WS there
                    cd ../3_Corridor_Skeleton_Made
                    save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS');
                    % Navigate to Video_Processed, remove the skeleton fields
                    % in WS, and then save WS there
                    cd ../2_Video_Processed
                    WS = rmfield(WS, {'LocOfSmallestGapInOddFlip','Skeleton_x','Skeleton_y'});
                    save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS');
                    disp(['Saved: ', pwd,'\', ...
                          WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'])
                else
                    % Navigate to Video_Processed and save WS once the video's processed
                    cd ../2_Video_Processed
                    save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS');
                    disp(['Saved: ', pwd,'\', ...
                          WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'])
                end
            end

        case MakeSkeleton
            % Loop through all experiments loaded in
            for vidCounter = 1:NumGCExp

                % Load in WS for each experiment
                WS = load(WS_names{vidCounter});
                % To get the correct WS, we must go in one layer
                WS = WS.WS;

                % Move to the appropriate folder for this experiment
                cd([LocalGCDirectoryPath,'Data\',...
                    WS.directoryName,'\Experiment_', num2str(WS.ExpRepNum)])

                % Navigate to Corridor_Skeleton_Made folder
                cd 3_Corridor_Skeleton_Made

                % Run CorridorSkeletonFinder
                WS = CorridorSkeletonFinderFunc(WS);

                % Save WS in Skeleton_Made
                save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS');
                disp(['Saved: ', pwd,'\', ...
                      WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'])
            end

        case RunNeuralNet
            % Loop through all experiments loaded in
            for vidCounter = 1:NumGCExp
                % Load in WS for each experiment
                WS = load(WS_names{vidCounter});
                % To get the correct WS, we must go in one layer
                WS = WS.WS;

                % Move to the appropriate folder for this experiment
                cd([LocalGCDirectoryPath,'Data\',...
                    WS.directoryName,'\Experiment_', num2str(WS.ExpRepNum)])

                % Navigate to Neural_Net_Analyzed folder
                cd 4_Neural_Net_Analyzed

                % Run all functions between finding skeleton and applying net
                WS = SkeletonToCorrMasksFunc(WS);
                WS = SkeletonToGapCentersFunc(WS);
                WS = SkeletonToCompMasksFunc(WS);
                WS = FlyCorrIdentifierFunc(WS);
                WS = FlyCompIdentifierFunc(WS);
                WS = FinalStatsToFlyStructFunc(WS);
                WS = FlyActivityFilterFunc(WS);
                WS = MetaDataAdderFunc(WS);
                WS = FlipBinnerFunc(WS);
                WS = ZapStationaryFlipsFunc(WS);
                WS = FindCrossEventsFunc(WS);
                WS = CrossingEventFrameFinderFunc(WS);
                WS = PredictOnOffGlassFunc(WS);

                % Save WS in Neural_Net_Analyzed
                save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS');
                disp(['Saved: ', pwd,'\', ...
                      WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'])
            end

        case RunFinalStep
            % Loop through all experiments loaded in
            for vidCounter = 1:NumGCExp
                % Load in WS for each experiment
                WS = load(WS_names{vidCounter});
                % To get the correct WS, we must go in one layer
                WS = WS.WS;

                % Move to the appropriate folder for this experiment
                cd([LocalGCDirectoryPath,'Data\',...
                    WS.directoryName,'\Experiment_', num2str(WS.ExpRepNum)])

                % Navigate to Fully_Analyzed
                cd 5_Fully_Analyzed

                % Run the last few functions to clean FBFS and give
                % classifications given the predicted probabilities
                WS = FinalFlipRemoverFunc(WS);
                WS = ClassifyOnOffFunc(WS);

                % Save WS in Fully_Analyzed and inform user
                save([WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'],'WS');
                disp(['Saved: ', pwd,'\', ...
                      WS.directoryName,'_Exp',num2str(WS.ExpRepNum),'_WS_struct.mat'])

                % Navigate to and save FlipBinnedFlyStruct in Fly_Structure
                cd ../Fly_Structure
                FlipBinnedFlyStruct = WS.FlipBinnedFlyStruct;
                save([WS.directoryName,'_FlipBinnedFlyStruct_Exp',num2str(WS.ExpRepNum),'.mat'],...
                      'FlipBinnedFlyStruct');
            end
            
        case AnalysisDone
            fprintf('All of the files have already been fully analyzed.\n')
            % Return to the folder that contains the analysis code and exit function
            cd([LocalGCDirectoryPath,'Matlab_Analysis_Scripts\FlipperRobotCode'])
            return

    end
    
    % Now progress to the next step by updating WS_names to Next_WS_names
    WS_names = Next_WS_names;
    
toc

    % If finished running the final step, inform user and return out of this function
    if strcmp(AnalysisStep, RunFinalStep)
        fprintf('All of the files have been fully analyzed.\n')
        % Return to the folder that contains the analysis code and exit function
        cd([LocalGCDirectoryPath,'Matlab_Analysis_Scripts\FlipperRobotCode'])
        return
    end

    % Give user choice to stop here since everything has now been saved,
    % but only do this if the user indicated that they wanted checkpoints
    if ~strcmp(availCheckpoints,'No')
        SaveCheckpointResponse = SaveCheckpointPromptFunc(contOpt, stopOpt);
    end

end

% Return to the folder that contains the analysis code
cd([LocalGCDirectoryPath,'Matlab_Analysis_Scripts\FlipperRobotCode'])

% Inform user that the analysis progress has been saved in command window
% once they choose to stop the analysis partially completed
fprintf('Your analysis progress has been saved.\n');

end