%% 1. Find datapaths to be processed
% Get data paths from the server by specifying stimulus parameter, cell
% type, indicators, etc.

% your cell type
% cellType = 'T4T5';
cellType = 'LC15';

% your indicator
sensor = 'GC6f';

% your name
surgeon = 'Braedyn';

% your stimulus
% stim = 'sinewave_20deg_3sEpoch_4sInt_4reps';
% stim = 'sinewave_30deg_fullCont_3sEpoch_4sInt_4reps';
% stim = 'sinewave_60deg_fullCont_3sEpoch_4sInt_4reps';
stim = 'bar40_80_160_walking';

% the eye you did the experiment
% this can remain empty unless you are ambidextrous...
flyEye = '';

% get path to the data
dataPath = GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
% Get flyIDs
flyIDs = GetFlyIdFromDatabase(dataPath);
%% 2. Go through each experiment and run AnalyzeFicTracDataBlockAverage
expResults = cell(length(dataPath),1);
% interleaveEpoch = 13;
% interleaveDur = 240;
% numEpochs = 14;
interleaveEpoch = 2;
interleaveDur = 300; % in block averaged behavior, this buffer becomes 30 frames at 14Hz
buffer = interleaveDur/10;
numEpochs = 18;
forceFictrac = 1;
redoProjections = 1;
% for ff = 53:length(dataPath)
% goodPaths = [12:15, 27:30 , 34:37 , 41:52, 54:57,63:66];
% goodPaths = [12:15, 27:30 , 34:37 , 41:56,63:66];
goodPaths = [12:15, 27:30 , 34:37 , 41:52, 54:57, 63:66]; % 53 is bad on elizabeth's computer so go back to original
% analyze and threshold
for ff = goodPaths
    disp([num2str(ff) '/' num2str(length(dataPath))]);
    if isfile([dataPath{ff} '\stimulusData\fictracData.dat'])
        disp(['Fictrac data found for recording ' num2str(ff) '/' num2str(length(dataPath))]);
        % doing block averaging
        [expResults{ff}.fictracAnalysis, expResults{ff}.allWalkingSpeeds, ...
            expResults{ff}.allTurningSpeeds, expResults{ff}.allForwardSpeeds, ... 
            expResults{ff}.processedData] = AnalyzeFicTracDataBlockAverage(dataPath{ff}, ...
            interleaveEpoch, interleaveDur, forceFictrac, redoProjections);
        expResults{ff}.flyID = flyIDs(ff);
    else % fictrac was not recording
        disp(['FICTRAC NOT ON FOR RECORDING ' num2str(ff) '/' num2str(length(dataPath)) ])
    end
end
%% 3. Go through all the analyzed data, and find the activity threshold for each unique fly across experiments
%% Right now it is filtering active/inactive based on inst. speed. Change to turning/forward walking here or later
% we have turning speed, forward speed, and inst. walking speed
% for fly = unique(flyIDs)
% goodPaths = [12:15, 27:30 , 34:37 , 41:56,63:66];
goodDataPaths = dataPath(goodPaths);
goodFlyIDs = flyIDs(goodPaths);
% for fly = unique(flyIDs(goodPaths))
% for fly =13393
for fly = unique(flyIDs(goodPaths))
    dataPathInds = find(flyIDs==fly);
    thisFlyAllTurningSpeeds = [];
    thisFlyAllWalkingSpeeds = [];
    thisFlyAllForwardSpeeds = [];
    for ff = dataPathInds
        if isempty(expResults{ff})
        % fictrac not recorded
            continue
        else
            thisFlyAllTurningSpeeds = [thisFlyAllTurningSpeeds; expResults{ff}.allTurningSpeeds]; % deg/s
            thisFlyAllWalkingSpeeds = [thisFlyAllWalkingSpeeds; expResults{ff}.allWalkingSpeeds]; % mm/s
            thisFlyAllForwardSpeeds = [thisFlyAllForwardSpeeds; expResults{ff}.allForwardSpeeds]; % mm/s
        end  
    end
    sphereRadius = 4.5; %mm
    mmToDeg = 180/pi/sphereRadius;
    % otsu threshold 
    % graythresh requires normalized values and returns a normalized
    % threshold so it needs to be mulitplied by max() to be applied
    % walkingThreshold = graythresh(rescale(thisFlyAllWalkingSpeeds))*max(thisFlyAllWalkingSpeeds);
    walkingThreshold = 1.5;
    % will this be too noisy? Take into account backwards walking?
    turningThreshold = graythresh(rescale(abs(thisFlyAllTurningSpeeds)))*max(abs(thisFlyAllTurningSpeeds));
    forwardThreshold = graythresh(rescale(abs(thisFlyAllForwardSpeeds)))*max(abs(thisFlyAllForwardSpeeds));
    % plot data 14 Hz after 10 block averaging
    timeX = linspace(0,length(thisFlyAllWalkingSpeeds)/14,length(thisFlyAllWalkingSpeeds));
    figure; hold on; set(gcf,'Position',[442, 462, 1085, 491]); t = tiledlayout(3,3); title(t,['Fly: ', num2str(fly)]);
    nexttile(1,[1,2]); hold on; plot(timeX,thisFlyAllWalkingSpeeds,'HandleVisibility','off'); yline(walkingThreshold, 'LineWidth',1.5,'Color','red'); xlabel('Time (s)');ylabel('Inst. Speed (mm/s)');
    nexttile(4,[1,2]); hold on; plot(1:length(thisFlyAllForwardSpeeds),thisFlyAllForwardSpeeds,'HandleVisibility','off'); yline(forwardThreshold, 'LineWidth',1.5,'Color','red'); xlabel('Time (s)');ylabel('Forward Walking Speed (mm/s)');
    nexttile(7,[1,2]); hold on; plot(1:length(thisFlyAllTurningSpeeds),thisFlyAllTurningSpeeds,'HandleVisibility','off'); yline(turningThreshold, 'LineWidth',1.5,'Color','red'); xlabel('Time (s)');ylabel('Turning Speed (deg/s)');
    % histograms
    nexttile(3); hold on; histogram(thisFlyAllWalkingSpeeds,'HandleVisibility','off'); set(gca,'YScale','log'); 
    xline(walkingThreshold, 'LineWidth',1.5,'Color','red'), xlabel('Inst. Speed (mm/s)'); ylabel('Count'); legend('Otsu Threshold');
    nexttile(6); hold on; histogram(thisFlyAllForwardSpeeds,'HandleVisibility','off'); set(gca,'YScale','log'); 
    xline(forwardThreshold, 'LineWidth',1.5,'Color','red'), xline(-forwardThreshold, 'LineWidth',1.5,'Color','red'); xlabel('Forward Walking Speed (mm/s)'); ylabel('Count'); legend('Otsu Threshold');%xlim([-mmToDeg*5, mmToDeg*5]);
    nexttile(9); hold on; histogram(thisFlyAllTurningSpeeds,'HandleVisibility','off'); set(gca,'YScale','log'); 
    xline(turningThreshold, 'LineWidth',1.5,'Color','red'), xline(-turningThreshold, 'LineWidth',1.5,'Color','red'); xlabel('Turning Speed (deg/s)'); ylabel('Count'); legend('Otsu Threshold');%xlim([-250, 250]);
    
    % save the thresholds in the dataPath folder
    if ~isempty(walkingThreshold) && ~isempty(turningThreshold) && ~isempty(forwardThreshold)
        for ff = dataPathInds
            if isempty(expResults{ff})
                continue
            end
            save([dataPath{ff} '\stimulusData\walkingThreshold.mat'], 'walkingThreshold');
            save([dataPath{ff} '\stimulusData\turningThreshold.mat'], 'turningThreshold');
            save([dataPath{ff} '\stimulusData\forwardThreshold.mat'], 'forwardThreshold');
            %% seperate active and quiescent trials
            % save a file called activeQuiescent with numEpochs x 1 cell
            % categorizing each trial based on inst. walking speed
            % other methods are done afterwards in plot...
            activeQuiescent = cell(numEpochs+interleaveEpoch,1);
            for ep = interleaveEpoch + 1:interleaveEpoch + numEpochs
                % start a boolean cell of active/quiescent trials
                numTrials = length(expResults{ff}.fictracAnalysis{ep}.trialWalkingSpeed);
                activeQuiescent{ep} = zeros(numTrials,1);
                for tr = 1:numTrials
                    numPoints = length(expResults{ff}.fictracAnalysis{ep}.trialWalkingSpeed{tr});
                    activeCounts = sum(expResults{ff}.fictracAnalysis{ep}.trialWalkingSpeed{tr}(buffer:numPoints-buffer) > walkingThreshold);
                    epochCounts = length(expResults{ff}.fictracAnalysis{ep}.trialWalkingSpeed{tr}(buffer:numPoints-buffer));
                    if activeCounts/epochCounts >= 0.25
                        % active trial b/c at least 25% activity throughout
                        activeQuiescent{ep}(tr) = 1;
                    end
                end
            end
    % save the finished activeQuiescent cell to each dataPath to be loaded
            save([dataPath{ff} '\stimulusData\activeQuiescent.mat'], 'activeQuiescent');
            save([dataPath{ff} '\stimulusData\thisFlyAllTurningSpeeds.mat'], 'thisFlyAllTurningSpeeds');
            save([dataPath{ff} '\stimulusData\thisFlyAllWalkingSpeeds.mat'], 'thisFlyAllWalkingSpeeds');
            save([dataPath{ff} '\stimulusData\thisFlyAllForwardSpeeds.mat'], 'thisFlyAllForwardSpeeds');
        end
    end
end
%% 4. Plot Behavior and Neuron response
roiExtractionFile = 'ManualRoiExtraction';

% use correlation based thresholding
% this is not a cleanly written function -- consider improving...        
%roiSelectionFile = 'selectROIbyProbeCorrelationGeneric'; 
roiSelectionFile = ''; % or don't use selection at alls


% analysisFiles = 'PlotTimeTraces'; % just plot trial averaged time traces
analysisFiles = 'PlotActiveInactiveTimeTracesElizabeth';

% dataPath([12:15, 27:30 , 34:37 , 41:52, 54:57, 63:66])
% goodPaths = [12:15, 27:30 , 34:37 , 41:56,63:66]; % this one or above? 53
% is bad on elizabeth's computer so its the one above
% goodPaths = [12:15, 27:30 , 34:37 , 41:52, 54:57, 63:66]; June 20, 2025
args = {'analysisFile',analysisFiles,...
        'dataPath',dataPath(goodPaths),...
        'roiExtractionFile',roiExtractionFile,...
        'roiSelectionFile',roiSelectionFile,...
        'forceRois',0,... % if set to 1, redo ROI extraction
        'individual',1,...
        'backgroundSubtractMovie', 1,...
        'backgroundSubtractByRoi', 0,...
        'calcDFOverFByRoi',1,...
        'epochsForSelectivity',{'dummy'},...
        'combOpp',0,...
        'epochsForIdentificationForFly',1,...
        'stimulusResponseAlignment',0,...
        'noTrueInterleave',0,...
        'perRoiDfOverFCalcFunction','CalculatedDeltaFOverFByROI',...
        'numPresentationForSelection',0,...
        'legacyFitting',0};
a1 = RunAnalysis(args{:});
%% 5. FWHM fitting for quiescent data
%% Quiescent
clrs = get(gca,'colororder');

mean_resp_mat = cell2mat(a1.analysis{1, 1}.quiescentRespPlot.');
num_flies = size(a1.analysis{1, 1}.indFly,2);

T = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Braedyn\LC15 Experiments\bar40_80_160_walking.txt');
% T = readcell('C:\Users\clare\Documents\Braedyn\psycho5\paramfiles\Braedyn\LC15 Experiments\bar40_80_160_walking.txt');
epoch_names = T(3,3:end);
epoch_durations = T(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);

bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;

x_ax_timings = a1.analysis{1, 1}.timeX ;
x_ax_timings_truncated = zeros(size(x_ax_timings,1), num_epochs);
for epoch = 1:num_epochs
    timing_in_epoch = x_ax_timings;
    idx_not_in_epoch = x_ax_timings > epoch_durations(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated(:,epoch) = timing_in_epoch;
end

quiescent_flies_avg_resp_fwhm = zeros(num_flies, num_epochs);
quiescent_flies_fwhm_lower_time = zeros(num_flies, num_epochs);
quiescent_flies_fwhm_upper_time = zeros(num_flies, num_epochs);

%calculate avg resp
for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        if length(a1.analysis{1, 1}.indFly{1, fly_num}.p9_quiescentTrialResps{epoch, 1}) < 1
            continue
        end
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p11_quiescentResps{epoch, 1};
    end
    

    fwhm_ind_fly = zeros(1, num_epochs);
    fwhm_lower_time_ind_fly = zeros(1, num_epochs);
    fwhm_upper_time_ind_fly = zeros(1, num_epochs);

    figure()
    hold on
    first = 1;
    for bar40_idx = bar40_idxs
        if (~all(ind_mean_resp_mat(:,bar40_idx) == 0) && (sum(isnan(ind_mean_resp_mat(:,bar40_idx))) < 0.1*length(ind_mean_resp_mat(:,bar40_idx))))
            curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
            curr_epoch_idxs = ~isnan(curr_epoch_idxs);
        
            curr_x_ax = x_ax_timings;
            curr_x_ax = curr_x_ax(curr_epoch_idxs);
        
            curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
            curr_y_ax = curr_y_ax(curr_epoch_idxs);
            
            if first == 1
                PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar40_idx,:));
            
                prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
                bounds_msec = inputdlg(prompt);
                bounds_msec = str2num(bounds_msec{1});
            
                lower_bound_msec = bounds_msec(1);
                upper_bound_msec = bounds_msec(2);
                first = 0;
            end
    %         lower_bound_msec = 4500;
    %         upper_bound_msec = 8500;
        
            %convert input bounds (in msec) to sampling rate that directly
            %correspond to resp mat idxs
            lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
            upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
        
            resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
            resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
        
            halfmax = (min(resp_window)+max(resp_window))/2;
            halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
            halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
            if halfmax_idx_sampled_lower_minus1 >= 1
                while isnan(resp_window(halfmax_idx_sampled_lower_minus1))
                    halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower_minus1 - 1;
                end
            end
            halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
            halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                while isnan(resp_window(halfmax_idx_sampled_upperplus1))
                    halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upperplus1 + 1;
                end
            end
    
          %interpolate to find two times of halfmax
    
            %avoid OOB (lower)
            if halfmax_idx_sampled_lower_minus1 >= 1
                halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_lower_minus1 = 1;
                halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
            end
          
            %avoid OOB (upper)
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
                halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
            end
    
    
            %insert these interpolated HM points into the vectors
            resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
            
            resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
                                            resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
                                            resp_window(halfmax_idx_sampled_upperplus1:end)];
    
            halfmax_idx_lower = halfmax_idx_sampled_lower;
            halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
    
            %integrate
            area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
            delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
            mean_value_resp = area/delta_time_sec;
        
            fwhm_ind_fly(bar40_idx) = mean_value_resp;
            fwhm_lower_time_ind_fly(bar40_idx) = halfmax_time_lower;
            fwhm_upper_time_ind_fly(bar40_idx) = halfmax_time_upper;
        else
            fwhm_ind_fly(bar40_idx) = NaN;
            fwhm_lower_time_ind_fly(bar40_idx) = NaN;
            fwhm_upper_time_ind_fly(bar40_idx) = NaN;
        end
    end
    hold off

    figure()
    hold on
    first = 1;
    for bar80_idx = bar80_idxs
        if (~all(ind_mean_resp_mat(:,bar80_idx) == 0) && (sum(isnan(ind_mean_resp_mat(:,bar80_idx))) < 0.1*length(ind_mean_resp_mat(:,bar80_idx))))
            curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
            curr_epoch_idxs = ~isnan(curr_epoch_idxs);
        
            curr_x_ax = x_ax_timings;
            curr_x_ax = curr_x_ax(curr_epoch_idxs);
        
            curr_y_ax = ind_mean_resp_mat(:,bar80_idx);
            curr_y_ax = curr_y_ax(curr_epoch_idxs);
            if first == 1
                PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar80_idx-6,:));
            
                prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
                bounds_msec = inputdlg(prompt);
                bounds_msec = str2num(bounds_msec{1});
            
                lower_bound_msec = bounds_msec(1);
                upper_bound_msec = bounds_msec(2);
                first = first +1;
            end
            
    %         lower_bound_msec = 2500;
    %         upper_bound_msec = 5500;
            %convert input bounds (in msec) to sampling rate that directly
            %correspond to resp mat idxs
            lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
            upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
        
            resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
            resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
        
            halfmax = (min(resp_window)+max(resp_window))/2;
            halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
            halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
            if halfmax_idx_sampled_lower_minus1 >= 1
                while isnan(resp_window(halfmax_idx_sampled_lower_minus1))
                    halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower_minus1 - 1;
                end
            end
            halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
            halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                while isnan(resp_window(halfmax_idx_sampled_upperplus1))
                    halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upperplus1 + 1;
                end
            end
    
    
          %interpolate to find two times of halfmax
    
            %avoid OOB (lower)
            if halfmax_idx_sampled_lower_minus1 >= 1
                halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_lower_minus1 = 1;
                halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
            end
          
            %avoid OOB (upper)
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
                halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
            end
    
    
            %insert these interpolated HM points into the vectors
            resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
            
            resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
                                            resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
                                            resp_window(halfmax_idx_sampled_upperplus1:end)];
    
            halfmax_idx_lower = halfmax_idx_sampled_lower;
            halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
    
            %integrate
            area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
            delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
            mean_value_resp = area/delta_time_sec;
        
            fwhm_ind_fly(bar80_idx) = mean_value_resp;
            fwhm_lower_time_ind_fly(bar80_idx) = halfmax_time_lower;
            fwhm_upper_time_ind_fly(bar80_idx) = halfmax_time_upper;
        else
            fwhm_ind_fly(bar80_idx) = NaN;
            fwhm_lower_time_ind_fly(bar80_idx) = NaN;
            fwhm_upper_time_ind_fly(bar80_idx) = NaN;
        end
    end
    hold off

    figure()
    hold on
    first = 1;
    for bar160_idx = bar160_idxs
        if (~all(ind_mean_resp_mat(:,bar160_idx) == 0) && (sum(isnan(ind_mean_resp_mat(:,bar160_idx))) < 0.1*length(ind_mean_resp_mat(:,bar160_idx))))
            curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
            curr_epoch_idxs = ~isnan(curr_epoch_idxs);
        
            curr_x_ax = x_ax_timings;
            curr_x_ax = curr_x_ax(curr_epoch_idxs);
        
            curr_y_ax = ind_mean_resp_mat(:,bar160_idx);
            curr_y_ax = curr_y_ax(curr_epoch_idxs);
            
            if first == 1
                PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar160_idx-12,:));
            
                prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
                bounds_msec = inputdlg(prompt);
                bounds_msec = str2num(bounds_msec{1});
            
                lower_bound_msec = bounds_msec(1);
                upper_bound_msec = bounds_msec(2);
                first = 0;
            end
    
    %         lower_bound_msec = 1500;
    %         upper_bound_msec = 4000;
            %convert input bounds (in msec) to sampling rate that directly
            %correspond to resp mat idxs
            lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
            upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
        
            resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
            resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
        
            halfmax = (min(resp_window)+max(resp_window))/2;
            halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
            halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
            if halfmax_idx_sampled_lower_minus1 >= 1
                while isnan(resp_window(halfmax_idx_sampled_lower_minus1))
                    halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower_minus1 - 1;
                end
            end
            halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
            halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                while isnan(resp_window(halfmax_idx_sampled_upperplus1))
                    halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upperplus1 + 1;
                end
            end
    
          %interpolate to find two times of halfmax
    
            %avoid OOB (lower)
            if halfmax_idx_sampled_lower_minus1 >= 1
                halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_lower_minus1 = 1;
                halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
            end
          
            %avoid OOB (upper)
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
                halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
            end
    
            
    
            %insert these interpolated HM points into the vectors
            resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
            
            resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
                                            resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
                                            resp_window(halfmax_idx_sampled_upperplus1:end)];
    
            halfmax_idx_lower = halfmax_idx_sampled_lower;
            halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
    
            %integrate
            area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
            delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
            mean_value_resp = area/delta_time_sec;
        
            fwhm_ind_fly(bar160_idx) = mean_value_resp;
            fwhm_lower_time_ind_fly(bar160_idx) = halfmax_time_lower;
            fwhm_upper_time_ind_fly(bar160_idx) = halfmax_time_upper;
        else
            fwhm_ind_fly(bar160_idx) = NaN;
            fwhm_lower_time_ind_fly(bar160_idx) = NaN;
            fwhm_upper_time_ind_fly(bar160_idx) = NaN;
        end
    end
    hold off

    quiescent_flies_avg_resp_fwhm(fly_num,:) = fwhm_ind_fly;
    quiescent_flies_fwhm_lower_time(fly_num,:) = fwhm_lower_time_ind_fly;
    quiescent_flies_fwhm_upper_time(fly_num,:) = fwhm_upper_time_ind_fly;
end

%% plotting quiescent
ratio = false;

figure();


v_ratio_x_ax = [0.0625 0.25 0.5 1 2 4];
bar40_fwhms_ratios = quiescent_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_bckg0 = bar40_fwhms_ratios(:, 1);

if ratio
    bar40_fwhms_ratios = bar40_fwhms_ratios./bar40_fwhms_bckg0;
    %bar40_fwhms_ratios = bar40_fwhms_ratios(:,2:end);
end

bar40_fwhms_ratios_mean = nanmean(bar40_fwhms_ratios, 1);

bar40_fwhms_ratios_sem = nanstd(bar40_fwhms_ratios, 1)/sqrt(num_flies);


bar80_fwhms_ratios = quiescent_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_bckg0 = bar80_fwhms_ratios(:, 1);

if ratio
    bar80_fwhms_ratios = bar80_fwhms_ratios./bar80_fwhms_bckg0;
    %bar80_fwhms_ratios = bar80_fwhms_ratios(:,2:end);
end

bar80_fwhms_ratios_mean = nanmean(bar80_fwhms_ratios, 1);

bar80_fwhms_ratios_sem = nanstd(bar80_fwhms_ratios, 1)/sqrt(num_flies);

bar160_fwhms_ratios = quiescent_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_bckg0 = bar160_fwhms_ratios(:, 1);

if ratio
    bar160_fwhms_ratios = bar160_fwhms_ratios./bar160_fwhms_bckg0;
    %bar160_fwhms_ratios = bar160_fwhms_ratios(:,2:end);
end

bar160_fwhms_ratios_mean = nanmean(bar160_fwhms_ratios, 1);

bar160_fwhms_ratios_sem = nanstd(bar160_fwhms_ratios, 1)/sqrt(num_flies);
hAx=axes;
hAx.XScale='log';   
hold all
errorbar(v_ratio_x_ax, (bar40_fwhms_ratios_mean), (bar40_fwhms_ratios_sem), 'Color', [0 0 0]);
errorbar(v_ratio_x_ax, (bar80_fwhms_ratios_mean), (bar80_fwhms_ratios_sem), 'Color', [0.75 0 0]);
errorbar(v_ratio_x_ax, (bar160_fwhms_ratios_mean), (bar160_fwhms_ratios_sem), 'Color', [0 0 0.75]);
xlim([1/24 4.5])
ylim([0 6])
xticks(v_ratio_x_ax(2:end))
xlabel('v_{background}/v_{bar}');
ylabel(['Mean Response - ', num2str(num_flies), ' Flies'])
title('Quiescent')
if ratio
    ylabel('Mean response (normalized)')
end

legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})

figure();
hold on
bar40_bckgvels = ([2.5 10 20 40 80 160]);
bar80_bckgvels = ([2.5 20 40 80 160 320]);
bar160_bckgvels = ([2.5 40 80 160 320 640]);


bar40_fwhms = quiescent_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_mean = nanmean(bar40_fwhms, 1);
bar40_fwhms_sem = nanstd(bar40_fwhms, 1)/sqrt(num_flies);

bar80_fwhms = quiescent_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_mean = nanmean(bar80_fwhms, 1);
bar80_fwhms_sem = nanstd(bar80_fwhms, 1)/sqrt(num_flies);

bar160_fwhms = quiescent_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_mean = nanmean(bar160_fwhms, 1);
bar160_fwhms_sem = nanstd(bar160_fwhms, 1)/sqrt(num_flies);

if ratio
    errorbar(bar40_bckgvels, bar40_fwhms_ratios_mean, bar40_fwhms_ratios_sem, 'Color', [0 0 0]);
    errorbar(bar80_bckgvels, bar80_fwhms_ratios_mean, bar80_fwhms_ratios_sem, 'Color', [0.75 0 0]);
    errorbar(bar160_bckgvels, bar160_fwhms_ratios_mean, bar160_fwhms_ratios_sem, 'Color', [0 0 0.75]);
else
    errorbar(bar40_bckgvels, bar40_fwhms_mean, bar40_fwhms_sem, 'Color', [0 0 0]);
    errorbar(bar80_bckgvels, bar80_fwhms_mean, bar80_fwhms_sem, 'Color', [0.75 0 0]);
    errorbar(bar160_bckgvels, bar160_fwhms_mean, bar160_fwhms_sem, 'Color', [0 0 0.75]);
end


xlabel('v_{bckg}');
ylabel(['Mean Response - ', num2str(num_flies), ' Flies'])
if ratio
    ylabel('Mean response (normalized)')
end

legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})
set(gca, 'XScale', 'log')
ylim([0 6])
title('Quiescent')
hold off
xticks([10,20,40,80,160,320,640])

%% Active
clrs = get(gca,'colororder');

mean_resp_mat = cell2mat(a1.analysis{1, 1}.activeRespPlot.');
num_flies = size(a1.analysis{1, 1}.indFly,2);

T = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Braedyn\LC15 Experiments\bar40_80_160_walking.txt');
% T = readcell('C:\Users\clare\Documents\Braedyn\psycho5\paramfiles\Braedyn\LC15 Experiments\bar40_80_160_walking.txt');
epoch_names = T(3,3:end);
epoch_durations = T(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);

bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;

x_ax_timings = a1.analysis{1, 1}.timeX ;
x_ax_timings_truncated = zeros(size(x_ax_timings,1), num_epochs);
for epoch = 1:num_epochs
    timing_in_epoch = x_ax_timings;
    idx_not_in_epoch = x_ax_timings > epoch_durations(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated(:,epoch) = timing_in_epoch;
end

active_flies_avg_resp_fwhm = zeros(num_flies, num_epochs);
active_flies_fwhm_lower_time = zeros(num_flies, num_epochs);
active_flies_fwhm_upper_time = zeros(num_flies, num_epochs);

%calculate avg resp
for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        if length(a1.analysis{1, 1}.indFly{1, fly_num}.p8_activeTrialResps{epoch, 1}) < 1
            continue
        end
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p10_activeResps{epoch, 1};
    end
    

    fwhm_ind_fly = zeros(1, num_epochs);
    fwhm_lower_time_ind_fly = zeros(1, num_epochs);
    fwhm_upper_time_ind_fly = zeros(1, num_epochs);

    figure()
    hold on
    first = 1;
    for bar40_idx = bar40_idxs
        if (~all(ind_mean_resp_mat(:,bar40_idx) == 0) && (sum(isnan(ind_mean_resp_mat(:,bar40_idx))) < 0.1*length(ind_mean_resp_mat(:,bar40_idx))))
            curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
            curr_epoch_idxs = ~isnan(curr_epoch_idxs);
        
            curr_x_ax = x_ax_timings;
            curr_x_ax = curr_x_ax(curr_epoch_idxs);
        
            curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
            curr_y_ax = curr_y_ax(curr_epoch_idxs);
            
            if first == 1
                PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar40_idx,:));
            
                prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
                bounds_msec = inputdlg(prompt);
                bounds_msec = str2num(bounds_msec{1});
            
                lower_bound_msec = bounds_msec(1);
                upper_bound_msec = bounds_msec(2);
                first = 0;
            end
    %         lower_bound_msec = 4500;
    %         upper_bound_msec = 8500;
        
            %convert input bounds (in msec) to sampling rate that directly
            %correspond to resp mat idxs
            lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
            upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
        
            resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
            resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
        
            halfmax = (min(resp_window)+max(resp_window))/2;
            halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
            halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
            if halfmax_idx_sampled_lower_minus1 >= 1
                while isnan(resp_window(halfmax_idx_sampled_lower_minus1))
                    halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower_minus1 - 1;
                end
            end
            halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
            halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                while isnan(resp_window(halfmax_idx_sampled_upperplus1))
                    halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upperplus1 + 1;
                end
            end
    
          %interpolate to find two times of halfmax
    
            %avoid OOB (lower)
            if halfmax_idx_sampled_lower_minus1 >= 1
                halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_lower_minus1 = 1;
                halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
            end
          
            %avoid OOB (upper)
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
                halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
            end
    
    
            %insert these interpolated HM points into the vectors
            resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
            
            resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
                                            resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
                                            resp_window(halfmax_idx_sampled_upperplus1:end)];
    
            halfmax_idx_lower = halfmax_idx_sampled_lower;
            halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
    
            %integrate
            area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
            delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
            mean_value_resp = area/delta_time_sec;
        
            fwhm_ind_fly(bar40_idx) = mean_value_resp;
            fwhm_lower_time_ind_fly(bar40_idx) = halfmax_time_lower;
            fwhm_upper_time_ind_fly(bar40_idx) = halfmax_time_upper;
        else
            fwhm_ind_fly(bar40_idx) = NaN;
            fwhm_lower_time_ind_fly(bar40_idx) = NaN;
            fwhm_upper_time_ind_fly(bar40_idx) = NaN;
        end
    end
    hold off

    figure()
    hold on
    first = 1;
    for bar80_idx = bar80_idxs
        if (~all(ind_mean_resp_mat(:,bar80_idx) == 0) && (sum(isnan(ind_mean_resp_mat(:,bar80_idx))) < 0.1*length(ind_mean_resp_mat(:,bar80_idx))))
            curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
            curr_epoch_idxs = ~isnan(curr_epoch_idxs);
        
            curr_x_ax = x_ax_timings;
            curr_x_ax = curr_x_ax(curr_epoch_idxs);
        
            curr_y_ax = ind_mean_resp_mat(:,bar80_idx);
            curr_y_ax = curr_y_ax(curr_epoch_idxs);
            if first == 1
                PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar80_idx-6,:));
            
                prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
                bounds_msec = inputdlg(prompt);
                bounds_msec = str2num(bounds_msec{1});
            
                lower_bound_msec = bounds_msec(1);
                upper_bound_msec = bounds_msec(2);
                first = first +1;
            end
            
    %         lower_bound_msec = 2500;
    %         upper_bound_msec = 5500;
            %convert input bounds (in msec) to sampling rate that directly
            %correspond to resp mat idxs
            lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
            upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
        
            resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
            resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
        
            halfmax = (min(resp_window)+max(resp_window))/2;
            halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
            halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
            if halfmax_idx_sampled_lower_minus1 >= 1
                while isnan(resp_window(halfmax_idx_sampled_lower_minus1))
                    halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower_minus1 - 1;
                end
            end
            halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
            halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                while isnan(resp_window(halfmax_idx_sampled_upperplus1))
                    halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upperplus1 + 1;
                end
            end
    
    
          %interpolate to find two times of halfmax
    
            %avoid OOB (lower)
            if halfmax_idx_sampled_lower_minus1 >= 1
                halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_lower_minus1 = 1;
                halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
            end
          
            %avoid OOB (upper)
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
                halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
            end
    
    
            %insert these interpolated HM points into the vectors
            resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
            
            resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
                                            resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
                                            resp_window(halfmax_idx_sampled_upperplus1:end)];
    
            halfmax_idx_lower = halfmax_idx_sampled_lower;
            halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
    
            %integrate
            area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
            delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
            mean_value_resp = area/delta_time_sec;
        
            fwhm_ind_fly(bar80_idx) = mean_value_resp;
            fwhm_lower_time_ind_fly(bar80_idx) = halfmax_time_lower;
            fwhm_upper_time_ind_fly(bar80_idx) = halfmax_time_upper;
        else
            fwhm_ind_fly(bar80_idx) = NaN;
            fwhm_lower_time_ind_fly(bar80_idx) = NaN;
            fwhm_upper_time_ind_fly(bar80_idx) = NaN;
        end
    end
    hold off

    figure()
    hold on
    first = 1;
    for bar160_idx = bar160_idxs
        if (~all(ind_mean_resp_mat(:,bar160_idx) == 0) && (sum(isnan(ind_mean_resp_mat(:,bar160_idx))) < 0.1*length(ind_mean_resp_mat(:,bar160_idx))))
            curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
            curr_epoch_idxs = ~isnan(curr_epoch_idxs);
        
            curr_x_ax = x_ax_timings;
            curr_x_ax = curr_x_ax(curr_epoch_idxs);
        
            curr_y_ax = ind_mean_resp_mat(:,bar160_idx);
            curr_y_ax = curr_y_ax(curr_epoch_idxs);
            
            if first == 1
                PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar160_idx-12,:));
            
                prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
                bounds_msec = inputdlg(prompt);
                bounds_msec = str2num(bounds_msec{1});
            
                lower_bound_msec = bounds_msec(1);
                upper_bound_msec = bounds_msec(2);
                first = 0;
            end
    
    %         lower_bound_msec = 1500;
    %         upper_bound_msec = 4000;
            %convert input bounds (in msec) to sampling rate that directly
            %correspond to resp mat idxs
            lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
            upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
        
            resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
            resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
        
            halfmax = (min(resp_window)+max(resp_window))/2;
            halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
            halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
            if halfmax_idx_sampled_lower_minus1 >= 1
                while isnan(resp_window(halfmax_idx_sampled_lower_minus1))
                    halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower_minus1 - 1;
                end
            end
            halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
            halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                while isnan(resp_window(halfmax_idx_sampled_upperplus1))
                    halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upperplus1 + 1;
                end
            end
    
          %interpolate to find two times of halfmax
    
            %avoid OOB (lower)
            if halfmax_idx_sampled_lower_minus1 >= 1
                halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_lower_minus1 = 1;
                halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
            end
          
            %avoid OOB (upper)
            if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
                halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
                                         [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
                                         halfmax);
            else
                halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
                halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
            end
    
            
    
            %insert these interpolated HM points into the vectors
            resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
                                            resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
            
            resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
                                            resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
                                            resp_window(halfmax_idx_sampled_upperplus1:end)];
    
            halfmax_idx_lower = halfmax_idx_sampled_lower;
            halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
    
            %integrate
            area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
            delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
            mean_value_resp = area/delta_time_sec;
        
            fwhm_ind_fly(bar160_idx) = mean_value_resp;
            fwhm_lower_time_ind_fly(bar160_idx) = halfmax_time_lower;
            fwhm_upper_time_ind_fly(bar160_idx) = halfmax_time_upper;
        else
            fwhm_ind_fly(bar160_idx) = NaN;
            fwhm_lower_time_ind_fly(bar160_idx) = NaN;
            fwhm_upper_time_ind_fly(bar160_idx) = NaN;
        end
    end
    hold off

    active_flies_avg_resp_fwhm(fly_num,:) = fwhm_ind_fly;
    active_flies_fwhm_lower_time(fly_num,:) = fwhm_lower_time_ind_fly;
    active_flies_fwhm_upper_time(fly_num,:) = fwhm_upper_time_ind_fly;
end

% %% Active
% clrs = get(gca,'colororder');
% 
% mean_resp_mat = cell2mat(a1.analysis{1, 1}.activeRespPlot.');
% num_flies = size(a1.analysis{1, 1}.indFly,2);
% 
% % T = readcell('C:\Users\clark\Documents\MATLAB\psycho5\paramfiles\Braedyn\LC15 Experiments\bar40_80_160_walking.txt');
% epoch_names = T(3,3:end);
% epoch_durations = T(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
% num_epochs  = size(epoch_names, 2);
% 
% bar40_idxs = 1:6;
% bar80_idxs = 7:12;
% bar160_idxs = 13:18;
% 
% x_ax_timings = a1.analysis{1, 1}.timeX ;
% x_ax_timings_truncated = zeros(size(x_ax_timings,1), num_epochs);
% for epoch = 1:num_epochs
%     timing_in_epoch = x_ax_timings;
%     idx_not_in_epoch = x_ax_timings > epoch_durations(epoch);
%     timing_in_epoch(idx_not_in_epoch) = NaN;
%     x_ax_timings_truncated(:,epoch) = timing_in_epoch;
% end
% 
% active_flies_avg_resp_fwhm = zeros(num_flies, num_epochs);
% active_flies_fwhm_lower_time = zeros(num_flies, num_epochs);
% active_flies_fwhm_upper_time = zeros(num_flies, num_epochs);
% 
% %calculate avg resp
% for fly_num = 1:num_flies
%     ind_mean_resp_mat = zeros(size(mean_resp_mat));
%     for epoch = 1:num_epochs
%         ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p10_activeResps{epoch, 1};
%     end
% 
% 
%     fwhm_ind_fly = zeros(1, num_epochs);
%     fwhm_lower_time_ind_fly = zeros(1, num_epochs);
%     fwhm_upper_time_ind_fly = zeros(1, num_epochs);
% 
%     figure()
%     hold on
%     first = 1;
%     for bar40_idx = bar40_idxs
%         curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
%         curr_epoch_idxs = ~isnan(curr_epoch_idxs);
% 
%         curr_x_ax = x_ax_timings;
%         curr_x_ax = curr_x_ax(curr_epoch_idxs);
% 
%         curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
%         curr_y_ax = curr_y_ax(curr_epoch_idxs);
%         if all(isnan(curr_y_ax))
%             fwhm_ind_fly(bar40_idx) = NaN;
%             fwhm_lower_time_ind_fly(bar40_idx) = NaN;
%             fwhm_upper_time_ind_fly(bar40_idx) = NaN;
%             continue
%         end
% 
%         if first
%             PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar40_idx,:));
% 
%             prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
%             bounds_msec = inputdlg(prompt);
%             bounds_msec = str2num(bounds_msec{1});
% 
%             lower_bound_msec = bounds_msec(1);
%             upper_bound_msec = bounds_msec(2);
%             first = 0;
%         end
% %         lower_bound_msec = 5000;
% %         upper_bound_msec = 7500;
% 
%         %convert input bounds (in msec) to sampling rate that directly
%         %correspond to resp mat idxs
% 
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
%         upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
% 
%         resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
%         resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
% 
%         halfmax = (min(resp_window)+max(resp_window))/2;
%         halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
%         halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
%         halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
%         halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
% 
%       %interpolate to find two times of halfmax
% 
%         %avoid OOB (lower)
%         if halfmax_idx_sampled_lower_minus1 >= 1
%             halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
%                                      [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
%                                      halfmax);
%         else
%             halfmax_idx_sampled_lower_minus1 = 1;
%             halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
%         end
% 
%         %avoid OOB (upper)
%         if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
%             halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
%                                      [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
%                                      halfmax);
%         else
%             halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
%             halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
%         end
% 
% 
%         %insert these interpolated HM points into the vectors
%         resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
%                                         resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
%                                         resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
% 
%         resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
%                                         resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
%                                         resp_window(halfmax_idx_sampled_upperplus1:end)];
% 
%         halfmax_idx_lower = halfmax_idx_sampled_lower;
%         halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
% 
%         %integrate
%         area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
%         delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
%         mean_value_resp = area/delta_time_sec;
% 
%         fwhm_ind_fly(bar40_idx) = mean_value_resp;
%         fwhm_lower_time_ind_fly(bar40_idx) = halfmax_time_lower;
%         fwhm_upper_time_ind_fly(bar40_idx) = halfmax_time_upper;
%     end
%     hold off
% 
%     figure()
%     hold on
%     first = 1;
%     for bar80_idx = bar80_idxs
%         curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
%         curr_epoch_idxs = ~isnan(curr_epoch_idxs);
% 
%         curr_x_ax = x_ax_timings;
%         curr_x_ax = curr_x_ax(curr_epoch_idxs);
% 
%         curr_y_ax = ind_mean_resp_mat(:,bar80_idx);
%         curr_y_ax = curr_y_ax(curr_epoch_idxs);
%         if all(isnan(curr_y_ax))
%             fwhm_ind_fly(bar80_idx) = NaN;
%             fwhm_lower_time_ind_fly(bar80_idx) = NaN;
%             fwhm_upper_time_ind_fly(bar80_idx) = NaN;
%             continue
%         end
% 
%         if first
%             PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar80_idx-6,:));
% 
%             prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
%             bounds_msec = inputdlg(prompt);
%             bounds_msec = str2num(bounds_msec{1});
% 
%             lower_bound_msec = bounds_msec(1);
%             upper_bound_msec = bounds_msec(2);
%             first = 0;
%         end
% 
% %         lower_bound_msec = 3000;
% %         upper_bound_msec = 4500;
%         %convert input bounds (in msec) to sampling rate that directly
%         %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
%         upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
% 
%         resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
%         resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
% 
%         halfmax = (min(resp_window)+max(resp_window))/2;
%         halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
%         halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
%         halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
%         halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
% 
%       %interpolate to find two times of halfmax
% 
%         %avoid OOB (lower)
%         if halfmax_idx_sampled_lower_minus1 >= 1
%             halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
%                                      [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
%                                      halfmax);
%         else
%             halfmax_idx_sampled_lower_minus1 = 1;
%             halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
%         end
% 
%         %avoid OOB (upper)
%         if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
%             halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
%                                      [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
%                                      halfmax);
%         else
%             halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
%             halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
%         end
% 
% 
%         %insert these interpolated HM points into the vectors
%         resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
%                                         resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
%                                         resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
% 
%         resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
%                                         resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
%                                         resp_window(halfmax_idx_sampled_upperplus1:end)];
% 
%         halfmax_idx_lower = halfmax_idx_sampled_lower;
%         halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
% 
%         %integrate
%         area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
%         delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
%         mean_value_resp = area/delta_time_sec;
% 
%         fwhm_ind_fly(bar80_idx) = mean_value_resp;
%         fwhm_lower_time_ind_fly(bar80_idx) = halfmax_time_lower;
%         fwhm_upper_time_ind_fly(bar80_idx) = halfmax_time_upper;
%     end
%     hold off
% 
%     figure()
%     hold on
%     first = 1;
%     for bar160_idx = bar160_idxs
%         curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
%         curr_epoch_idxs = ~isnan(curr_epoch_idxs);
% 
%         curr_x_ax = x_ax_timings;
%         curr_x_ax = curr_x_ax(curr_epoch_idxs);
% 
%         curr_y_ax = ind_mean_resp_mat(:,bar160_idx);
%         curr_y_ax = curr_y_ax(curr_epoch_idxs);
%         if all(isnan(curr_y_ax))
%             fwhm_ind_fly(bar160_idx) = NaN;
%             fwhm_lower_time_ind_fly(bar160_idx) = NaN;
%             fwhm_upper_time_ind_fly(bar160_idx) = NaN;
%             continue
%         end
%         if first
%             PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar160_idx-12,:));
% 
%             prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
%             bounds_msec = inputdlg(prompt);
%             bounds_msec = str2num(bounds_msec{1});
% 
%             lower_bound_msec = bounds_msec(1);
%             upper_bound_msec = bounds_msec(2);
%             first = 0;
%         end
% 
% %         lower_bound_msec = 2000;
% %         upper_bound_msec = 3000;
%         %convert input bounds (in msec) to sampling rate that directly
%         %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
%         upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
% 
%         resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
%         resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
% 
%         halfmax = (min(resp_window)+max(resp_window))/2;
%         halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
%         halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
%         halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
%         halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;
% 
%       %interpolate to find two times of halfmax
% 
%         %avoid OOB (lower)
%         if halfmax_idx_sampled_lower_minus1 >= 1
%             halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
%                                      [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
%                                      halfmax);
%         else
%             halfmax_idx_sampled_lower_minus1 = 1;
%             halfmax_time_lower = resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1);
%         end
% 
%         %avoid OOB (upper)
%         if halfmax_idx_sampled_upperplus1 <= size(resp_window,1)
%             halfmax_time_upper = interp1([resp_window(halfmax_idx_sampled_upper) resp_window(halfmax_idx_sampled_upperplus1)],...
%                                      [resp_window_curr_x_ax(halfmax_idx_sampled_upper) resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1)],...
%                                      halfmax);
%         else
%             halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper;
%             halfmax_time_upper = resp_window_curr_x_ax(halfmax_idx_sampled_upper);
%         end
% 
% 
% 
%         %insert these interpolated HM points into the vectors
%         resp_window_curr_x_ax_interp = [resp_window_curr_x_ax(1:halfmax_idx_sampled_lower_minus1); halfmax_time_lower;...
%                                         resp_window_curr_x_ax(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax_time_upper; ...
%                                         resp_window_curr_x_ax(halfmax_idx_sampled_upperplus1:end)];
% 
%         resp_window_interp = [resp_window(1:halfmax_idx_sampled_lower_minus1); halfmax;...
%                                         resp_window(halfmax_idx_sampled_lower:halfmax_idx_sampled_upper); halfmax; ...
%                                         resp_window(halfmax_idx_sampled_upperplus1:end)];
% 
%         halfmax_idx_lower = halfmax_idx_sampled_lower;
%         halfmax_idx_upper = halfmax_idx_sampled_upper + 2; 
% 
%         %integrate
%         area = trapz(resp_window_curr_x_ax_interp(halfmax_idx_lower:halfmax_idx_upper), resp_window_interp(halfmax_idx_lower:halfmax_idx_upper)); 
%         delta_time_sec = (halfmax_time_upper - halfmax_time_lower);
%         mean_value_resp = area/delta_time_sec;
% 
%         if ~isempty(mean_value_resp)
%             fwhm_ind_fly(bar160_idx) = mean_value_resp;
%             fwhm_lower_time_ind_fly(bar160_idx) = halfmax_time_lower;
%             fwhm_upper_time_ind_fly(bar160_idx) = halfmax_time_upper;
%         else
%             fwhm_ind_fly(bar160_idx) = nan;
%             fwhm_lower_time_ind_fly(bar160_idx) = nan;
%             fwhm_upper_time_ind_fly(bar160_idx) = nan;
%         end
%     end
%     hold off
% 
%     active_flies_avg_resp_fwhm(fly_num,:) = fwhm_ind_fly;
%     active_flies_fwhm_lower_time(fly_num,:) = fwhm_lower_time_ind_fly;
%     active_flies_fwhm_upper_time(fly_num,:) = fwhm_upper_time_ind_fly;
% end

%% plotting active
ratio = false;

figure();

v_ratio_x_ax = [0.0625 0.25 0.5 1 2 4];
bar40_fwhms_ratios = active_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_bckg0 = bar40_fwhms_ratios(:, 1);

if ratio
    bar40_fwhms_ratios = bar40_fwhms_ratios./bar40_fwhms_bckg0;
    %bar40_fwhms_ratios = bar40_fwhms_ratios(:,2:end);
end

bar40_fwhms_ratios_mean = nanmean(bar40_fwhms_ratios, 1);

bar40_fwhms_ratios_sem = nanstd(bar40_fwhms_ratios, 1)/sqrt(num_flies);


bar80_fwhms_ratios = active_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_bckg0 = bar80_fwhms_ratios(:, 1);

if ratio
    bar80_fwhms_ratios = bar80_fwhms_ratios./bar80_fwhms_bckg0;
    %bar80_fwhms_ratios = bar80_fwhms_ratios(:,2:end);
end

bar80_fwhms_ratios_mean = nanmean(bar80_fwhms_ratios, 1);

bar80_fwhms_ratios_sem = nanstd(bar80_fwhms_ratios, 1)/sqrt(num_flies);

bar160_fwhms_ratios = active_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_bckg0 = bar160_fwhms_ratios(:, 1);

if ratio
    bar160_fwhms_ratios = bar160_fwhms_ratios./bar160_fwhms_bckg0;
    %bar160_fwhms_ratios = bar160_fwhms_ratios(:,2:end);
end

bar160_fwhms_ratios_mean = nanmean(bar160_fwhms_ratios, 1);

bar160_fwhms_ratios_sem = nanstd(bar160_fwhms_ratios, 1)/sqrt(num_flies);
hAx=axes;
hAx.XScale='log';   
hold all
errorbar(v_ratio_x_ax, (bar40_fwhms_ratios_mean), (bar40_fwhms_ratios_sem), 'Color', [0 0 0]);
errorbar(v_ratio_x_ax, (bar80_fwhms_ratios_mean), (bar80_fwhms_ratios_sem), 'Color', [0.75 0 0]);
errorbar(v_ratio_x_ax, (bar160_fwhms_ratios_mean), (bar160_fwhms_ratios_sem), 'Color', [0 0 0.75]);
xlim([1/24 4.5])
ylim([0 6])
xticks(v_ratio_x_ax(2:end))
xlabel('v_{background}/v_{bar}');
ylabel(['Mean Response - ', num2str(num_flies), ' Flies'])
title('Active')
if ratio
    ylabel('Mean response (normalized)')
end

legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})

figure();
hold on
bar40_bckgvels = ([2.5 10 20 40 80 160]);
bar80_bckgvels = ([2.5 20 40 80 160 320]);
bar160_bckgvels = ([2.5 40 80 160 320 640]);


bar40_fwhms = active_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_mean = nanmean(bar40_fwhms, 1);
bar40_fwhms_sem = nanstd(bar40_fwhms, 1)/sqrt(num_flies);

bar80_fwhms = active_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_mean = nanmean(bar80_fwhms, 1);
bar80_fwhms_sem = nanstd(bar80_fwhms, 1)/sqrt(num_flies);

bar160_fwhms = active_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_mean = nanmean(bar160_fwhms, 1);
bar160_fwhms_sem = nanstd(bar160_fwhms, 1)/sqrt(num_flies);

if ratio
    errorbar(bar40_bckgvels, bar40_fwhms_ratios_mean, bar40_fwhms_ratios_sem, 'Color', [0 0 0]);
    errorbar(bar80_bckgvels, bar80_fwhms_ratios_mean, bar80_fwhms_ratios_sem, 'Color', [0.75 0 0]);
    errorbar(bar160_bckgvels, bar160_fwhms_ratios_mean, bar160_fwhms_ratios_sem, 'Color', [0 0 0.75]);
else
    errorbar(bar40_bckgvels, bar40_fwhms_mean, bar40_fwhms_sem, 'Color', [0 0 0]);
    errorbar(bar80_bckgvels, bar80_fwhms_mean, bar80_fwhms_sem, 'Color', [0.75 0 0]);
    errorbar(bar160_bckgvels, bar160_fwhms_mean, bar160_fwhms_sem, 'Color', [0 0 0.75]);
end


xlabel('v_{bckg}');
ylabel(['Mean Response - ', num2str(num_flies), ' Flies'])
if ratio
    ylabel('Mean response (normalized)')
end

legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})
set(gca, 'XScale', 'log')
ylim([0 6])
title('Active')
hold off
xticks([10,20,40,80,160,320,640])

%% bar plot active/quiescent difference
active_mean_resp = nanmean(active_flies_avg_resp_fwhm,1);
quiescent_mean_resp = nanmean(quiescent_flies_avg_resp_fwhm,1);
flyColorMap = prism(num_flies);
figure; hold on;
bp = bar(1:18,[quiescent_mean_resp;active_mean_resp]);
bp(1).FaceColor = [0.3,0.3,0.8];
bp(2).FaceColor = [0.8,0.3,0.3];
for ep = 1:18
    plot([bp(1).XEndPoints(ep),bp(2).XEndPoints(ep)],[quiescent_flies_avg_resp_fwhm(:,ep),active_flies_avg_resp_fwhm(:,ep)],'-k','HandleVisibility','off');
end
for fly = 1:num_flies
    scatter([bp(1).XEndPoints,bp(2).XEndPoints],[quiescent_flies_avg_resp_fwhm(fly,:),active_flies_avg_resp_fwhm(fly,:)],36, 'filled','HandleVisibility','off');
end
% star any significance
sig = ttest(quiescent_flies_avg_resp_fwhm,active_flies_avg_resp_fwhm);
for i = 1:18
    if sig(i) == 1
        plot(i,9, '*k', 'HandleVisibility','off')
    end
end
xticks(1:18)
xticklabels(epoch_names)
ylabel(['Mean Response - ', num2str(num_flies), ' Flies'])
legend({'Quiescent','Active'})
ylim([0,9.5])
%% scatter mean response and mean walk speed for each epoch within the time window
figure; hold on;
tiledlayout('flow'); 
for ep = 1:numEpochs
    nexttile; hold on;
    % title(figLeg{numIgnore+ep});
    for ff = 1:numFlies
        meanWalk = nanmean(a1.analysis{1, 1}.indFly{ff}.p1_fictrac.meanWalkingSpeeds{ep}(30:end-30));
        meanDeltaFQuies = nanmean(a1.analysis{1, 1}.indFly{ff}.p9_quiescentTrialResps{ep}(30:end-30));
        meanDeltaFAct = nanmean(a1.analysis{1, 1}.indFly{ff}.p8_activeTrialResps{ep}(30:end-30));
        scatter(meanWalk,meanDeltaFQuies,'blue');
        scatter(meanWalk,meanDeltaFAct,'red');
    end
    xlabel('Mean Walking Speed in FWHM window (mm/s)')
    ylabel('Max \DeltaF / F in FWHM window')
    ylim([-0.5,10])
    xlim([0,15])
end
