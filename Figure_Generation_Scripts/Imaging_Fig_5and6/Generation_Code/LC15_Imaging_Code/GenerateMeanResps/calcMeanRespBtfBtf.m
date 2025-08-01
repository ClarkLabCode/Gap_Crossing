% your cell type
cellType = 'LC15';

% your indicator
sensor = 'GC6f';

% your name
surgeon = 'Elizabeth';

% your stimulus
stim = 'bar40_80_160_flipped';

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
clrs = get(gca,'colororder');
T = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll_flipped\bar40_80_160_flipped.txt');
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

%all_flies_avg_resp_fwhm = zeros(num_flies, num_epochs);
%all_flies_fwhm_lower_time = zeros(num_flies, num_epochs);
%all_flies_fwhm_upper_time = zeros(num_flies, num_epochs);

%calculate avg resp
for fly_num = 12:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p8_averagedRois.snipMat{epoch, 1};
    end
    

    fwhm_ind_fly = zeros(1, num_epochs);
    fwhm_lower_time_ind_fly = zeros(1, num_epochs);
    fwhm_upper_time_ind_fly = zeros(1, num_epochs);

    figure()
    hold on
    
    for bar40_idx = bar40_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar40_idx,:));
    
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
    end
    hold off

    figure()
    hold on
    for bar80_idx = bar80_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar80_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar80_idx-6,:));
    
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
    end
    hold off

    figure()
    hold on
    for bar160_idx = bar160_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar160_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        PlotXvsY(curr_x_ax, curr_y_ax, 'color', clrs(bar160_idx-12,:));
    
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
    end
    hold off

    all_flies_avg_resp_fwhm(fly_num,:) = fwhm_ind_fly;
    all_flies_fwhm_lower_time(fly_num,:) = fwhm_lower_time_ind_fly;
    all_flies_fwhm_upper_time(fly_num,:) = fwhm_upper_time_ind_fly;
end

%%
ratio = false;

figure();


v_ratio_x_ax = [0.0625 0.25 0.5 1 2 4];
bar40_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_bckg0 = bar40_fwhms_ratios(:, 1);
    
if ratio
    bar40_fwhms_ratios = bar40_fwhms_ratios./bar40_fwhms_bckg0;
    %bar40_fwhms_ratios = bar40_fwhms_ratios(:,2:end);
end

bar40_fwhms_ratios_mean = nanmean(bar40_fwhms_ratios, 1);

bar40_fwhms_ratios_sem = nanstd(bar40_fwhms_ratios, 1)/sqrt(num_flies);


bar80_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_bckg0 = bar80_fwhms_ratios(:, 1);

if ratio
    bar80_fwhms_ratios = bar80_fwhms_ratios./bar80_fwhms_bckg0;
    %bar80_fwhms_ratios = bar80_fwhms_ratios(:,2:end);
end

bar80_fwhms_ratios_mean = nanmean(bar80_fwhms_ratios, 1);

bar80_fwhms_ratios_sem = nanstd(bar80_fwhms_ratios, 1)/sqrt(num_flies);

bar160_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar160_idxs);
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
ylim([0 4.5])
xticks(v_ratio_x_ax(2:end))
xlabel('v_{background}/v_{bar}');
ylabel('Mean response')
if ratio
    ylabel('Mean response (normalized)')
    ylim([0 1.5])
end

legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})

figure();
hold on
bar40_bckgvels = ([2.5 10 20 40 80 160]);
bar80_bckgvels = ([2.5 20 40 80 160 320]);
bar160_bckgvels = ([2.5 40 80 160 320 640]);


bar40_fwhms = all_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_mean = nanmean(bar40_fwhms, 1);
bar40_fwhms_sem = nanstd(bar40_fwhms, 1)/sqrt(num_flies);

bar80_fwhms = all_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_mean = nanmean(bar80_fwhms, 1);
bar80_fwhms_sem = nanstd(bar80_fwhms, 1)/sqrt(num_flies);

bar160_fwhms = all_flies_avg_resp_fwhm(:, bar160_idxs);
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
ylabel('Mean Response')
ylim([0 4.5])
if ratio
    ylabel('Mean response (normalized)')
    ylim([0 1.5])
end

legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})
set(gca, 'XScale', 'log')
hold off
xticks([10,20,40,80,160,320,640])
