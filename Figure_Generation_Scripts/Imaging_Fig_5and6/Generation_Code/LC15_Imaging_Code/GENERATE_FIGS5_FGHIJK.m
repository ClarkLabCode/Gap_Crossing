%% plotting params
time_trace_y_bounds = [-0.5 4.5];
time_trace_lower_t = -2000;
trace_t_scale_bar = [-2000 0];
trace_y_scale_bar = [0 2];
tuning_curve_y_bounds = [-0.25 4.5];
tuning_curve_yticks = [0 1 2 3 4];
bar_40_color = [0.8500 0.3250 0.0980];
bar_80_color = [0.4940 0.1840 0.5560];
bar_160_color = [0.3010 0.7450 0.9330];

%% S5F: bckg BTF bar FTB draw timeseries plots for rel. motion stims 


%get info for antidirectional stim
% your cell type
cellType = 'LC15';
% your indicator
sensor = 'GC6f';
% your name
surgeon = 'Elizabeth';
% your stimulus
stim = 'bar40_80_160_opposing';  
% the eye you did the experiment
% this can remain empty unless you are ambidextrous...
flyEye = '';
dataPath = GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
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

bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;
data_struct_anti = RunAnalysis(args{:});

%get epoch info from paramfile
num_flies_anti = size(data_struct_anti.analysis{1, 1}.indFly,2);
param_file_anti = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll\opposingMotion\bar40_80_160_opposing.txt');
epoch_names_anti = param_file_anti(3,3:end);
epoch_durations_anti = param_file_anti(7, 3:end); epoch_durations_anti = cell2mat(epoch_durations_anti); epoch_durations_anti = epoch_durations_anti./60*1000+1000;
num_epochs_anti  = size(epoch_names_anti, 2);

%get values for y axis
mean_resp_mat_anti = data_struct_anti.analysis{1, 1}.respMatPlot;
sem_resp_mat_anti = data_struct_anti.analysis{1, 1}.respMatSemPlot;

%make x-axis timings appropriate for bar speed
x_ax_timings_anti = data_struct_anti.analysis{1, 1}.timeX ;
x_ax_timings_truncated_anti = zeros(size(x_ax_timings_anti,1), num_epochs_anti);
for epoch = 1:num_epochs_anti
    timing_in_epoch = x_ax_timings_anti;
    idx_not_in_epoch = x_ax_timings_anti > epoch_durations_anti(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated_anti(:,epoch) = timing_in_epoch;
end


%3 bar speeds, 6 background speeds per
%each row is a bar speed [40 80 160]
%each col is bckg spd of ratio [0 .25 .5 1 2 4]
figure('Position', [1, 1, 600, 300]);
t = tiledlayout(3, 6,'TileSpacing','Tight');

for bar40_idx = bar40_idxs
    nexttile
    hold on

    xline(1000, "--k");
    ylim(time_trace_y_bounds);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar40_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar40_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar40_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', bar_40_color);

    xlim([time_trace_lower_t curr_x_ax(end)]);
    if bar40_idx == bar40_idxs(1)
        line(trace_t_scale_bar,[0,0], 'Color', [0 0 0])
        line([-2000 -2000],trace_y_scale_bar, 'Color', [0 0 0])
    end
    axis off
    hold off
end

for bar80_idx = bar80_idxs
    nexttile
    hold on

    xline(1000, "--k");
    ylim(time_trace_y_bounds);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar80_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar80_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar80_idx);
    curr_sem = curr_sem(curr_epoch_idxs);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    if bar80_idx == bar80_idxs(1)
        line(trace_t_scale_bar,[0,0], 'Color', [0 0 0])
        line([-2000 -2000],trace_y_scale_bar, 'Color', [0 0 0])
    end
    axis off
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', bar_80_color);
    hold off
end

for bar160_idx = bar160_idxs
    nexttile
    hold on

    xline(1000, "--k");
    ylim(time_trace_y_bounds);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar160_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar160_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar160_idx);
    curr_sem = curr_sem(curr_epoch_idxs);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    if bar160_idx == bar160_idxs(1)
        line(trace_t_scale_bar,[0,0], 'Color', [0 0 0])
        line([-2000 -2000],trace_y_scale_bar, 'Color', [0 0 0])
    end

    axis off
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', bar_160_color);
    hold off
end

%% S5G: plot mean resps wrt bg spd barFtb bgBtf
load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\opposing_relmot_means.mat");
num_flies = size(all_flies_avg_resp_fwhm,1);

%plot resp wrt backg spd (not ratio)
figure()
hold on
bar40_bckgvels = ([2.5 10 20 40 80 160]);
bar80_bckgvels = ([2.5 20 40 80 160 320]);
bar160_bckgvels = ([2.5 40 80 160 320 640]);


bar40_fwhms = all_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_mean = mean(bar40_fwhms, 1);
bar40_fwhms_sem = std(bar40_fwhms, 1)/sqrt(num_flies);

bar80_fwhms = all_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_mean = mean(bar80_fwhms, 1);
bar80_fwhms_sem = std(bar80_fwhms, 1)/sqrt(num_flies);

bar160_fwhms = all_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_mean = mean(bar160_fwhms, 1);
bar160_fwhms_sem = std(bar160_fwhms, 1)/sqrt(num_flies);

errorbar(bar40_bckgvels(2:end), (bar40_fwhms_mean(2:end)), (bar40_fwhms_sem(2:end)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
errorbar(bar80_bckgvels(2:end), (bar80_fwhms_mean(2:end)), (bar80_fwhms_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
errorbar(bar160_bckgvels(2:end), (bar160_fwhms_mean(2:end)), (bar160_fwhms_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330]);

errorbar(bar40_bckgvels(1), (bar40_fwhms_mean(1)), (bar40_fwhms_sem(1)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
errorbar(bar80_bckgvels(1), (bar80_fwhms_mean(1)), (bar80_fwhms_sem(1)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
errorbar(bar160_bckgvels(1), (bar160_fwhms_mean(1)), (bar160_fwhms_sem(1)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330]);


xlabel('v_{bckg}');
ylabel('Mean Response')
%legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
yticks(tuning_curve_yticks)
legend()
ylim(tuning_curve_y_bounds)
hold off
xticks([10,20,40,80,160,320,640])


%% S5H: plot wrt ratio of spds barFtb bgBtf

load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\opposing_relmot_means.mat");

figure();

v_ratio_x_ax = [0.0625 0.25 0.5 1 2 4];

bar40_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_bckg0 = bar40_fwhms_ratios(:, 1);
bar40_fwhms_ratios_mean = mean(bar40_fwhms_ratios, 1);
bar40_fwhms_ratios_sem = std(bar40_fwhms_ratios, 1)/sqrt(num_flies);

bar80_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_bckg0 = bar80_fwhms_ratios(:, 1);
bar80_fwhms_ratios_mean = mean(bar80_fwhms_ratios, 1);
bar80_fwhms_ratios_sem = std(bar80_fwhms_ratios, 1)/sqrt(num_flies);

bar160_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_bckg0 = bar160_fwhms_ratios(:, 1);
bar160_fwhms_ratios_mean = mean(bar160_fwhms_ratios, 1);
bar160_fwhms_ratios_sem = std(bar160_fwhms_ratios, 1)/sqrt(num_flies);

hAx=axes;
hAx.XScale='log';   
hold on

errorbar(v_ratio_x_ax(2:end), (bar40_fwhms_ratios_mean(2:end)), (bar40_fwhms_ratios_sem(2:end)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(2:end), (bar80_fwhms_ratios_mean(2:end)), (bar80_fwhms_ratios_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(2:end), (bar160_fwhms_ratios_mean(2:end)), (bar160_fwhms_ratios_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330], 'MarkerSize', 6);

errorbar(v_ratio_x_ax(1), (bar40_fwhms_ratios_mean(1)), (bar40_fwhms_ratios_sem(1)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(1), (bar80_fwhms_ratios_mean(1)), (bar80_fwhms_ratios_sem(1)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(1), (bar160_fwhms_ratios_mean(1)), (bar160_fwhms_ratios_sem(1)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330], 'MarkerSize', 6);

xlim([1/24 4.5])
ylim(tuning_curve_y_bounds)
set(gca,'XMinorTick','off')
yticks(tuning_curve_yticks)
xticks(v_ratio_x_ax(2:end))

xlabel('v_{background}/v_{bar}');
ylabel('Mean response')
legend({'v_{bar} = 40°/s', 'v_{bar} = 80°/s', 'v_{bar} = 160°/s'})
hold off


%% S5i: bckg BTF bar BTF draw timeseries plots for rel. motion stims 


%get info for antidirectional stim
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
dataPath = GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
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

bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;
data_struct_anti = RunAnalysis(args{:});

%get epoch info from paramfile
num_flies_anti = size(data_struct_anti.analysis{1, 1}.indFly,2);
param_file_anti = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll_flipped\bar40_80_160_flipped.txt');
epoch_names_anti = param_file_anti(3,3:end);
epoch_durations_anti = param_file_anti(7, 3:end); epoch_durations_anti = cell2mat(epoch_durations_anti); epoch_durations_anti = epoch_durations_anti./60*1000+1000;
num_epochs_anti  = size(epoch_names_anti, 2);

%get values for y axis
mean_resp_mat_anti = data_struct_anti.analysis{1, 1}.respMatPlot;
sem_resp_mat_anti = data_struct_anti.analysis{1, 1}.respMatSemPlot;

%make x-axis timings appropriate for bar speed
x_ax_timings_anti = data_struct_anti.analysis{1, 1}.timeX ;
x_ax_timings_truncated_anti = zeros(size(x_ax_timings_anti,1), num_epochs_anti);
for epoch = 1:num_epochs_anti
    timing_in_epoch = x_ax_timings_anti;
    idx_not_in_epoch = x_ax_timings_anti > epoch_durations_anti(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated_anti(:,epoch) = timing_in_epoch;
end


%3 bar speeds, 6 background speeds per
%each row is a bar speed [40 80 160]
%each col is bckg spd of ratio [0 .25 .5 1 2 4]
figure('Position', [1, 1, 600, 300]);
t = tiledlayout(3, 6,'TileSpacing','Tight');

for bar40_idx = bar40_idxs
    nexttile
    hold on

    xline(1000, "--k");
    ylim(time_trace_y_bounds);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar40_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar40_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar40_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', bar_40_color);

    xlim([time_trace_lower_t curr_x_ax(end)]);
    if bar40_idx == bar40_idxs(1)
        line(trace_t_scale_bar,[0,0], 'Color', [0 0 0])
        line([-2000 -2000],trace_y_scale_bar, 'Color', [0 0 0])
    end
    axis off
    hold off
end

for bar80_idx = bar80_idxs
    nexttile
    hold on

    xline(1000, "--k");
    ylim(time_trace_y_bounds);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar80_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar80_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar80_idx);
    curr_sem = curr_sem(curr_epoch_idxs);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    if bar80_idx == bar80_idxs(1)
        line(trace_t_scale_bar,[0,0], 'Color', [0 0 0])
        line([-2000 -2000],trace_y_scale_bar, 'Color', [0 0 0])
    end
    axis off
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', bar_80_color);
    hold off
end

for bar160_idx = bar160_idxs
    nexttile
    hold on

    xline(1000, "--k");
    ylim(time_trace_y_bounds);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar160_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar160_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar160_idx);
    curr_sem = curr_sem(curr_epoch_idxs);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    if bar160_idx == bar160_idxs(1)
        line(trace_t_scale_bar,[0,0], 'Color', [0 0 0])
        line([-2000 -2000],trace_y_scale_bar, 'Color', [0 0 0])
    end

    axis off
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', bar_160_color);
    hold off
end

%% S5j: plot mean resps wrt bg spd barBtf bgBtf

load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\flipped_relmot_means.mat");
num_flies = size(all_flies_avg_resp_fwhm,1);
%plot resp wrt backg spd (not ratio)
figure()
hold on
bar40_bckgvels = ([2.5 10 20 40 80 160]);
bar80_bckgvels = ([2.5 20 40 80 160 320]);
bar160_bckgvels = ([2.5 40 80 160 320 640]);


bar40_fwhms = all_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_mean = mean(bar40_fwhms, 1);
bar40_fwhms_sem = std(bar40_fwhms, 1)/sqrt(num_flies);

bar80_fwhms = all_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_mean = mean(bar80_fwhms, 1);
bar80_fwhms_sem = std(bar80_fwhms, 1)/sqrt(num_flies);

bar160_fwhms = all_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_mean = mean(bar160_fwhms, 1);
bar160_fwhms_sem = std(bar160_fwhms, 1)/sqrt(num_flies);

errorbar(bar40_bckgvels(2:end), (bar40_fwhms_mean(2:end)), (bar40_fwhms_sem(2:end)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
errorbar(bar80_bckgvels(2:end), (bar80_fwhms_mean(2:end)), (bar80_fwhms_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
errorbar(bar160_bckgvels(2:end), (bar160_fwhms_mean(2:end)), (bar160_fwhms_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330]);

errorbar(bar40_bckgvels(1), (bar40_fwhms_mean(1)), (bar40_fwhms_sem(1)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
errorbar(bar80_bckgvels(1), (bar80_fwhms_mean(1)), (bar80_fwhms_sem(1)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
errorbar(bar160_bckgvels(1), (bar160_fwhms_mean(1)), (bar160_fwhms_sem(1)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330]);


xlabel('v_{bckg}');
ylabel('Mean Response')
%legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
yticks(tuning_curve_yticks)
legend()
ylim(tuning_curve_y_bounds)
hold off
xticks([10,20,40,80,160,320,640])

%% S6K: plot mean resps wrt ratio of spds barBtf bgBtf

load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\flipped_relmot_means.mat");

figure();

v_ratio_x_ax = [0.0625 0.25 0.5 1 2 4];

bar40_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar40_idxs);
bar40_fwhms_bckg0 = bar40_fwhms_ratios(:, 1);
bar40_fwhms_ratios_mean = mean(bar40_fwhms_ratios, 1);
bar40_fwhms_ratios_sem = std(bar40_fwhms_ratios, 1)/sqrt(num_flies);

bar80_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar80_idxs);
bar80_fwhms_bckg0 = bar80_fwhms_ratios(:, 1);
bar80_fwhms_ratios_mean = mean(bar80_fwhms_ratios, 1);
bar80_fwhms_ratios_sem = std(bar80_fwhms_ratios, 1)/sqrt(num_flies);

bar160_fwhms_ratios = all_flies_avg_resp_fwhm(:, bar160_idxs);
bar160_fwhms_bckg0 = bar160_fwhms_ratios(:, 1);
bar160_fwhms_ratios_mean = mean(bar160_fwhms_ratios, 1);
bar160_fwhms_ratios_sem = std(bar160_fwhms_ratios, 1)/sqrt(num_flies);

hAx=axes;
hAx.XScale='log';   
hold on

errorbar(v_ratio_x_ax(2:end), (bar40_fwhms_ratios_mean(2:end)), (bar40_fwhms_ratios_sem(2:end)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(2:end), (bar80_fwhms_ratios_mean(2:end)), (bar80_fwhms_ratios_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(2:end), (bar160_fwhms_ratios_mean(2:end)), (bar160_fwhms_ratios_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330], 'MarkerSize', 6);

errorbar(v_ratio_x_ax(1), (bar40_fwhms_ratios_mean(1)), (bar40_fwhms_ratios_sem(1)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[0.8500 0.3250 0.0980], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(1), (bar80_fwhms_ratios_mean(1)), (bar80_fwhms_ratios_sem(1)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560], 'MarkerSize', 6);
errorbar(v_ratio_x_ax(1), (bar160_fwhms_ratios_mean(1)), (bar160_fwhms_ratios_sem(1)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[0.3010 0.7450 0.9330], 'MarkerSize', 6);

xlim([1/24 4.5])
ylim(tuning_curve_y_bounds)
set(gca,'XMinorTick','off')
yticks(tuning_curve_yticks)
xticks(v_ratio_x_ax(2:end))

xlabel('v_{background}/v_{bar}');
ylabel('Mean response')
legend({'v_{bar} = 40°/s', 'v_{bar} = 80°/s', 'v_{bar} = 160°/s'})
hold off




