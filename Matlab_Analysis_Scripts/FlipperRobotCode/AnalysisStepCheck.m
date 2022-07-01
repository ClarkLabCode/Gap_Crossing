function [AnalysisStep, Next_WS_names, SkeletonMade] = AnalysisStepCheck(WS_names)

% Create variables to track analysis step
ProcessVideo = 'Process Video';
MakeSkeleton = 'Make Skeleton';
RunNeuralNet = 'Run Neural Net';
RunFinalStep = 'Run Final Steps';
AnalysisDone = 'None';

% Check if WS comes from '1_Params_Entered' folder as proxy for
% checking if all parameters have been entered in the analysis
AnalysisStageOfExpVec = contains(WS_names, '1_Params_Entered');

% Check that all of the loaded files are in the same stage of parameter entry
if ~all(AnalysisStageOfExpVec == AnalysisStageOfExpVec(1))
    error('The loaded files are not at the same stage of analysis.')
end

% Once we know all files loaded are at the same stage of parameter entry,
% check whether or not they have had their parameters entered
if AnalysisStageOfExpVec(1)
    ParamsEntered = true;
else
    ParamsEntered = false;
end

% Check if WS comes from '2_Video_Processed' folder as proxy for
% checking if video processing has been completed in the analysis
AnalysisStageOfExpVec = contains(WS_names, '2_Video_Processed');

% Check that all of the loaded files are in the same stage of video processing
if ~all(AnalysisStageOfExpVec == AnalysisStageOfExpVec(1))
    error('The loaded files are not at the same stage of analysis.')
end

% Once we know all files loaded are at the same stage of video processing,
% check whether or not they have had their associated videos processed
if AnalysisStageOfExpVec(1)
    VidProcessed = true;
    % Have to also set previous analysis steps to true
    ParamsEntered = true;
else
    VidProcessed = false;
end

% Check if WS comes from '3_Corridor_Skeleton_Made' folder as proxy for
% checking if skeleton making has been completed in the analysis
AnalysisStageOfExpVec = contains(WS_names, '3_Corridor_Skeleton_Made');

% Check that all of the loaded files are in the same stage of skeleton making
if ~all(AnalysisStageOfExpVec == AnalysisStageOfExpVec(1))
    error('The loaded files are not at the same stage of analysis.')
end

% Once we know all files loaded are at the same stage of skeleton making,
% check whether or not they have had their skeletons made
if AnalysisStageOfExpVec(1)
    SkeletonMade = true;
    % Check if video processing has also been done by loading in one WS and
    % checking for finalStats field which is generated from video processing
    tempWS = load(['..\..\',WS_names{1}]);
    tempWS = tempWS.WS;
    % If it has, set VidProcessed to be true too
    if isfield(tempWS, 'finalStats')
        VidProcessed = true;
    else
        VidProcessed = false;
    end
    % Free up the memory occupied by tempWS
    clear tempWS
    % Have to also set previous analysis steps to true
    ParamsEntered = true;
else
    SkeletonMade = false;
end

% Check if WS comes from '4_Neural_Net_Analyzed' folder as proxy for
% checking if neural net feeding has been completed in the analysis
AnalysisStageOfExpVec = contains(WS_names, '4_Neural_Net_Analyzed');

% Check that all of the loaded files are in the same stage of neural net analysis
if ~all(AnalysisStageOfExpVec == AnalysisStageOfExpVec(1))
    error('The loaded files are not at the same stage of analysis.')
end

% Once we know all files loaded are at the same stage of neural net analysis,
% check whether or not they have been analyzed by the neural net
if AnalysisStageOfExpVec(1)
    NeuralNetAnalyzed = true;
    % Have to also set previous analysis steps to true
    SkeletonMade = true;
    VidProcessed = true;
    ParamsEntered = true;
else
    NeuralNetAnalyzed = false;
end

% Check if WS comes from '5_Fully_Analyzed' folder as proxy for
% checking if the files have been fully analyzed
AnalysisStageOfExpVec = contains(WS_names, '5_Fully_Analyzed');

% Check that all of the loaded files are in the same stage of full analysis
if ~all(AnalysisStageOfExpVec == AnalysisStageOfExpVec(1))
    error('The loaded files are not at the same stage of analysis.')
end

% Once we know all files loaded are at the same stage of full analysis,
% check whether or not they have been fully analyzed
if AnalysisStageOfExpVec(1)
    FullyAnalyzed = true;
    % Have to also set previous analysis steps to true
    NeuralNetAnalyzed = true;
    SkeletonMade = true;
    VidProcessed = true;
    ParamsEntered = true;
else
    FullyAnalyzed = false;
end

% If the experiments have been fully analyzed, set everything to done 
% and exit function
if FullyAnalyzed
    AnalysisStep = AnalysisDone;
    Next_WS_names = WS_names;
    return
end

% If the files don't come from any of the expected subfolders, throw an
% error to the user with some direction
if ~(ParamsEntered||VidProcessed||SkeletonMade||NeuralNetAnalyzed)
    error(['The files selected are not from the expected folders. ', ...
           'Please select .mat files only from one of the following folders:',...
           '\n   1_Params_Entered',...
           '\n   2_Video_Processed',...
           '\n   3_Corridor_Skeleton_Made',...
           '\n   4_Neural_Net_Analyzed'],1);
end

% If neither the video has been processed nor the skeleton has been made,
% give the user the option to choose which step to do first
if (~VidProcessed && ~SkeletonMade)
    AnalysisStep = questdlg('Which analysis step would you like to do first?',...
                                'Analysis Order', ProcessVideo, MakeSkeleton,...
                                ProcessVideo);
    % Based on what the user selected as the next analysis step, fill
    % Next_WS_names with the correct entry so that the step after next
    % can be properly determined by this function using Next_WS_names
    if strcmp(AnalysisStep, ProcessVideo)
        % In Next_WS_names, replace the folder the file came from with the
        % folder the new file was written into
        Next_WS_names = strrep(WS_names, '1_Params_Entered', '2_Video_Processed');
    elseif strcmp(AnalysisStep, MakeSkeleton)
        % In Next_WS_names, replace the folder the file came from with the
        % folder the new file was written into
        Next_WS_names = strrep(WS_names, '1_Params_Entered', '3_Corridor_Skeleton_Made');
    end
% If files have already been through the neural net, move to rest of analysis
elseif NeuralNetAnalyzed
    AnalysisStep = RunFinalStep;
    % In Next_WS_names, replace the folder the file came from with the
    % folder the new file was written into
    Next_WS_names = strrep(WS_names, '4_Neural_Net_Analyzed', '5_Fully_Analyzed');
% If files have had both the video processed and skeleton made but not been
% run through the neural net, move to the neural net analysis
elseif (VidProcessed && SkeletonMade)
    AnalysisStep = RunNeuralNet;
    % In Next_WS_names, replace the folder the file came from with the
    % folder the new file was written into
    Next_WS_names = strrep(WS_names, '3_Corridor_Skeleton_Made', '4_Neural_Net_Analyzed');
% If only the video has been processed, then move to make skeletong
elseif VidProcessed
    AnalysisStep = MakeSkeleton;
    % In Next_WS_names, replace the folder the file came from with the
    % folder the new file was written into
    Next_WS_names = strrep(WS_names, '2_Video_Processed', '3_Corridor_Skeleton_Made');
% If only the skeleton has been made, then move to video processing
elseif SkeletonMade
    AnalysisStep = ProcessVideo;
    % Because this step is out of the usual order, we don't change
    % Next_WS_names because the algorithm will figure it out on its own
    % with the same WS_names
    % This happens because the algorithm saves over the previous WS file
    % that had the skeleton made with a new WS file that still has the
    % skeleton made but also has the video processed, so when the "same"
    % file is loaded again, it passes both the video processing check and
    % the skeleton make check
    Next_WS_names = WS_names;
end

end
