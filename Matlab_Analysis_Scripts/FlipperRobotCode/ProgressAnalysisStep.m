function WS_names = ProgressAnalysisStep(WS_names)

% Depending on what the previous analysis step was, determine subsequent
% step and change WS_names accordingly so that AnalysisStepCheck.m can
% identify the correct subsequent step
switch PrevAnalysisStep
    % If video was processed just now, check whether or not skeleton
    % was already made. If it was, no need to change WS_names because
    % algorithm will find that skeleton and video are done and move to
    % the neural net step next. If it wasn't, then need to change from
    % 1_Params_Entered folder to 2_Video_Processed folder
    case ProcessVideo
        if ~SkeletonMade
            WS_names = strrep(WS_names,'1_Params_Entered','2_Video_Processed');
        end
    % If the skeleton was just made, need to change WS_names from the
    % previous folder it was in (either 1_Params_Entered or
    % 2_Video_Processed) to 3_Corridor_Skeleton_Made folder
    case MakeSkeleton
        if VidProcessed
            WS_names = strrep(WS_names,'2_Video_Processed','3_Corridor_Skeleton_Made');
        else
            WS_names = strrep(WS_names,'1_Params_Entered','3_Corridor_Skeleton_Made');
        end
    % If the neural net was just run, need to change WS_names from
    % 3_Corridor_Skeleton_Made folder to 4_Neural_Net_Analyzed folder
    case RunNeuralNet
        WS_names = strrep(WS_names,'3_Corridor_Skeleton_Made','4_Neural_Net_Analyzed');
    % If the final steps were just run, need to change WS_names from
    % 4_Neural_Net_Analyzed folder to 5_Fully_Analyzed folder
    case RunFinalStep
        WS_names = strrep(WS_names,'4_Neural_Net_Analyzed','5_Fully_Analyzed');
        % Also set FullyAnalyzed to true
        FullyAnalyzed = true;

end

end
