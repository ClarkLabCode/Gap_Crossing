% your cell type
cellType = 'LC15';

% your indicator
sensor = 'GC6f';

% your name
surgeon = 'Elizabeth';

% your stimulus
stim = 'bar20_40_80_160_320_LR_RL';

% the eye you did the experiment
% this can remain empty unless you are ambidextrous...
flyEye = '';

% get path to the data
dataPath = GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
%dataPath = 'Y:\2p_microscope_data\w_+;LC15AD_GC6f;LC15DBD_+\bar40_80_160_noOcclusion_ternary\2024\03_13\13_20_06';
roiExtractionFile = 'WatershedRegionRestrictedRoiExtraction'; % use watershed + manual circling around
roiExtractionFile = 'ManualRoiExtraction';

% use correlation based thresholding
% this is not a cleanly written function -- consider improving...        
%roiSelectionFile = 'selectROIbyProbeCorrelationGeneric'; 
roiSelectionFile = ''; % or don't use selection at alls


analysisFiles = 'PlotTimeTraces'; % just plot trial averaged time traces


args = {'analysisFile',analysisFiles,...
        'dataPath',dataPath,...
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
mean_resp_mat = a1.analysis{1, 1}.respMatPlot;
num_flies = size(a1.analysis{1, 1}.indFly,2);

T = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll\objSweep\bar20_40_80_160_320_LR_RL.txt');
epoch_names = T(3,3:end);
epoch_durations = T(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);

lr_idxs = 1:5;
rl_idxs = lr_idxs+5;

all_flies_avg_resp_fwhm = zeros(num_flies, num_epochs);
all_flies_fwhm_lower_time = zeros(num_flies, num_epochs);
all_flies_fwhm_upper_time = zeros(num_flies, num_epochs);

x_ax_timings = a1.analysis{1, 1}.timeX ;
x_ax_timings_truncated = zeros(size(x_ax_timings,1), num_epochs);
for epoch = 1:num_epochs
    timing_in_epoch = x_ax_timings;
    idx_not_in_epoch = x_ax_timings > epoch_durations(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated(:,epoch) = timing_in_epoch;
end
acq_rate = diff(x_ax_timings);
acq_rate = acq_rate(1)/1000; %in ms
 
for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p8_averagedRois.snipMat{epoch, 1};
    end
    

    fwhm_ind_fly = zeros(1, num_epochs);
    fwhm_lower_time_ind_fly = zeros(1, num_epochs);
    fwhm_upper_time_ind_fly = zeros(1, num_epochs);
    
    for barlr_idx = lr_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,barlr_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,barlr_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        PlotXvsY(curr_x_ax, curr_y_ax);
    
        prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
        bounds_msec = inputdlg(prompt);
        bounds_msec = str2num(bounds_msec{1});
    
        lower_bound_msec = bounds_msec(1);
        upper_bound_msec = bounds_msec(2);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
    
        halfmax = (min(resp_window)+max(resp_window))/2;
        halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
        halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
        halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
        halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;

        

        %interpolate to find two times of halfmax
        halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                     [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                     halfmax);

        %avoid OOB
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
    
        fwhm_ind_fly(barlr_idx) = mean_value_resp;
        fwhm_lower_time_ind_fly(barlr_idx) = halfmax_time_lower;
        fwhm_upper_time_ind_fly(barlr_idx) = halfmax_time_upper;
    end

    for barrl_idx = rl_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,barrl_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,barrl_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        PlotXvsY(curr_x_ax, curr_y_ax);
    
        prompt = "Enter lower and upper bounds to scan for FWHM (msec): ";
        bounds_msec = inputdlg(prompt);
        bounds_msec = str2num(bounds_msec{1});
    
        lower_bound_msec = bounds_msec(1);
        upper_bound_msec = bounds_msec(2);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        resp_window_curr_x_ax = curr_x_ax(lower_bound_idx:upper_bound_idx);
    
        halfmax = (min(resp_window)+max(resp_window))/2;
        halfmax_idx_sampled_lower = find(resp_window > halfmax, 1, "first");
        halfmax_idx_sampled_lower_minus1 = halfmax_idx_sampled_lower-1;
        halfmax_idx_sampled_upper = find(resp_window > halfmax, 1, "last");
        halfmax_idx_sampled_upperplus1 = halfmax_idx_sampled_upper+1;

        

        %interpolate to find two times of halfmax
        halfmax_time_lower = interp1([resp_window(halfmax_idx_sampled_lower) resp_window(halfmax_idx_sampled_lower_minus1)],...
                                     [resp_window_curr_x_ax(halfmax_idx_sampled_lower) resp_window_curr_x_ax(halfmax_idx_sampled_lower_minus1)],...
                                     halfmax);

        %avoid OOB
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
    
        fwhm_ind_fly(barrl_idx) = mean_value_resp;
        fwhm_lower_time_ind_fly(barrl_idx) = halfmax_time_lower;
        fwhm_upper_time_ind_fly(barrl_idx) = halfmax_time_upper;
    end
    
 

   

    all_flies_avg_resp_fwhm(fly_num,:) = fwhm_ind_fly;
    all_flies_fwhm_lower_time(fly_num,:) = fwhm_lower_time_ind_fly;
    all_flies_fwhm_upper_time(fly_num,:) = fwhm_upper_time_ind_fly;
end

%% 
figure();
num_flies =6;
x_ax = [20 40 80 160 320];

fwhm_means = mean(all_flies_avg_resp_fwhm,1);
fwhm_sems = std(all_flies_avg_resp_fwhm,1)/sqrt(num_flies);
hold on
errorbar(x_ax, fwhm_means(lr_idxs), fwhm_sems(lr_idxs), 'Color', [0 0 0]);
errorbar(x_ax, fwhm_means(rl_idxs), fwhm_sems(rl_idxs), 'Color', [0.75 0 0]);

xlabel('v_{bar}');
ylabel('Mean Response')
legend({'v_{bar} Progressive', 'v_{bar} Regressive'})
set(gca, 'XScale', 'log')
xlim([10 640])
xticks(x_ax);
ylim([0 6])
hold off
%%
anova_mat = zeros(0, 2);
for cond = lr_idxs
    anova_mat = [anova_mat; [all_flies_avg_resp_fwhm(:,cond) all_flies_avg_resp_fwhm(:,cond)]];

end
anova2(anova_mat,6)