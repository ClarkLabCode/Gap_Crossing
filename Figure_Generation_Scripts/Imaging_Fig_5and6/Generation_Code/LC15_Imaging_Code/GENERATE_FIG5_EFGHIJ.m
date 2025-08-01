%%plotting params
time_trace_y_bounds = [-0.5 4.5];
time_trace_lower_t = -2000;
trace_t_scale_bar = [-2000 0];
trace_y_scale_bar = [0 2];
tuning_curve_y_bounds = [-0.25 4.5];
tuning_curve_yticks = [0 1 2 3 4];
%% FIG 5E&F draw plots for obj (bar) sweep
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
data_struct = RunAnalysis(args{:});

%get epoch info from paramfile
num_flies = size(data_struct.analysis{1, 1}.indFly,2);
param_file = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll\objSweep\bar20_40_80_160_320_LR_RL.txt');
epoch_names = param_file(3,3:end);
epoch_durations = param_file(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);

lr_idxs = 1:5;
%note that the rl idx = is lr+5, this will be used in the plot loop


%get values for y axis
mean_resp_mat = data_struct.analysis{1, 1}.respMatPlot;
sem_resp_mat = data_struct.analysis{1, 1}.respMatSemPlot;

%make x-axis timings appropriate for bar speed
x_ax_timings = data_struct.analysis{1, 1}.timeX ;
x_ax_timings_truncated = zeros(size(x_ax_timings,1), num_epochs);
for epoch = 1:num_epochs
    timing_in_epoch = x_ax_timings;
    idx_not_in_epoch = x_ax_timings > epoch_durations(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated(:,epoch) = timing_in_epoch;
end

%5 bar speeds, 2 directions
%each col is bar spd [20 40 80 160 320]
figure('Position', [1, 1, 600, 100]);
t = tiledlayout(1, 5,'TileSpacing','Tight');

for bar_idx = lr_idxs
    nexttile
    curr_epoch_idxs = x_ax_timings_truncated(:,bar_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat(:,bar_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat(:,bar_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    %plot bckg start vs bar start

    ylim(time_trace_y_bounds);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem);

    hold on
    %note that the rl bar idx is lr_idx + 5
    curr_epoch_idxs = x_ax_timings_truncated(:,bar_idx+5);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat(:,bar_idx+5);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat(:,bar_idx+5);
    curr_sem = curr_sem(curr_epoch_idxs);
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [1 0.5 0.7] );
    %legend()
    
    line(trace_t_scale_bar,[0 0], 'Color', [0 0 0])
    line([-2000,-2000],trace_y_scale_bar, 'Color', [0 0 0])
    %xline(0, "--", Color = [.7 .7 .7]);
    xline(1000, "--k");
    axis off
    hold off;
end

% plot avg resps
load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\obj_sweep_LR_RL_means.mat");

figure();
num_flies = size(all_flies_fwhm, 1);
x_ax = [20 40 80 160 320];
rl_idxs = 6:10;

fwhm_means = mean(all_flies_fwhm,1);
fwhm_sems = std(all_flies_fwhm,1)/sqrt(num_flies);
hold on

errorbar(x_ax, fwhm_means(rl_idxs), fwhm_sems(rl_idxs),'o-', 'Color', [1 0.5 0.7], 'CapSize',0, 'MarkerFaceColor', [1 0.5 0.7]);
errorbar(x_ax, fwhm_means(lr_idxs), fwhm_sems(lr_idxs),'o-' , 'Color', [0 0.4470 0.7410], 'CapSize',0, 'MarkerFaceColor', [0 0.4470 0.7410]);

xlabel('v_{bar}');
ylabel('Mean Response')
legend({'v_{bar} Progressive', 'v_{bar} Regressive'})
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
xlim([10 640])
xticks(x_ax);
yticks(tuning_curve_yticks)
%line([25,25],[0,1], 'Color', [0 0 0])
ylim(tuning_curve_y_bounds)
hold off

% t-test for bar spd only - paired bc on same 6 flies
p_vals_bar = zeros(1, size(lr_idxs,2));
for bar_spd_idx = lr_idxs
    [h,p,ci,stats] = ttest(all_flies_fwhm(:,bar_spd_idx), all_flies_fwhm(:,bar_spd_idx+5)); %rl is lr +5
    p_vals_bar(bar_spd_idx) = p;
end

rej_null = find(p_vals_bar < 0.05);
p_vals_bar


%% FIG 5G GENERATE BCKG ONLY TRACES: not starts at -500 ms and time scalebar is 500 ms

% MAKE PROGRESSIVE BCKG RESP TRACES
% Get data paths from the server by specifying stimulus parameter, cell
% type, indicators, etc.

% your cell type
cellType = 'LC15';

% your indicator
sensor = 'GC6f';

% your name
surgeon = 'Elizabeth';

% your stimulus
stim = 'bar40_80_160';

% the eye you did the experiment
% this can remain empty unless you are ambidextrous...
flyEye = '';

% get path to the data
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
a1 = RunAnalysis(args{:});


mean_resp_mat = a1.analysis{1, 1}.respMatPlot;
num_flies = a1.analysis{1, 1}.numFlies;

%get stimulus parameters from paramfile
T = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll\bar40_80_160.txt');
epoch_names = T(3,3:end);
epoch_durations = T(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);

bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;

%truncate timeseries to only the relevant times
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
% lower_bound_msec = -2000; %only background displayed from 0 to 1 sec
lower_bound_idx = 1;
upper_bound_msec = 1000;
    
%loop through flies, then loops through timeseries of response to each
%epoch. Takes only the first 1000 msec of stim presentation because that is
%the time in which only the bckg is moving (bar has not started yet)
bckg_mean_resp_traces= cell(num_flies,1);

for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p8_averagedRois.snipMat{epoch, 1};
    end
    

    traces_ind_fly = zeros(25, num_epochs);

    
    for bar40_idx = bar40_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        traces_ind_fly(:, bar40_idx) = resp_window;


    end

    
    %BAR 80
    for bar80_idx = bar80_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar80_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        traces_ind_fly(:, bar80_idx) = resp_window;


    end


    
    %BAR 160
    for bar160_idx = bar160_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar160_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        traces_ind_fly(:, bar160_idx) = resp_window;

    end
    
    bckg_mean_resp_traces{fly_num} = traces_ind_fly;
end

bckg_spds = [0 10 20 40 80 160 320 640];
num_spds = length(bckg_spds);
mat_bckg_spds = repmat([0 10 20 40 80 160 0 20 40 80 160 320 0 40 80 160 320 640], num_flies,1);

mean_resp_by_bckg_spd_traces = cell(num_spds,1);
means_across_flies_resp_bckg_spd_prog_traces = zeros(25, num_spds);
sems_across_flies_resp_bckg_spd_prog_traces = zeros(25, num_spds);
for spd_idx = 1:num_spds
    curr_spd_resp_idxs = mat_bckg_spds == bckg_spds(spd_idx);
    curr_spd_resps = [];
    for fly = 1:num_flies
        curr_spd_resps = [curr_spd_resps, bckg_mean_resp_traces{fly}(:,curr_spd_resp_idxs(fly,:))];
    end
    mean_resp_by_bckg_spd_traces(spd_idx) = {curr_spd_resps};

    means_across_flies_resp_bckg_spd_prog_traces(:,spd_idx) = mean(curr_spd_resps,2);
    sems_across_flies_resp_bckg_spd_prog_traces(:,spd_idx) = std(curr_spd_resps,0,2)/sqrt(num_flies);
end

% BCKG REGRESSIVE
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

% get path to the data
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
a1 = RunAnalysis(args{:});


mean_resp_mat = a1.analysis{1, 1}.respMatPlot;
num_flies = a1.analysis{1, 1}.numFlies;

%get stimulus parameters from paramfile
T = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll\opposingMotion\bar40_80_160_opposing.txt');
epoch_names = T(3,3:end);
epoch_durations = T(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);

bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;

%truncate timeseries to only the relevant times
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
lower_bound_msec = -2000; %only background displayed from 0 to 1 sec
upper_bound_msec = 1000;
    
%loop through flies, then loops through timeseries of response to each
%epoch. Takes only the first 1000 msec of stim presentation because that is
%the time in which only the bckg is moving (bar has not started yet)
bckg_mean_resp_traces= cell(num_flies,1);

for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p8_averagedRois.snipMat{epoch, 1};
    end
    

    traces_ind_fly = zeros(25, num_epochs);

    
    for bar40_idx = bar40_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        traces_ind_fly(:, bar40_idx) = resp_window;


    end

    
    %BAR 80
    for bar80_idx = bar80_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar80_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        traces_ind_fly(:, bar80_idx) = resp_window;


    end


    
    %BAR 160
    for bar160_idx = bar160_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar160_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
%         lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
        traces_ind_fly(:, bar160_idx) = resp_window;

    end
    
    bckg_mean_resp_traces{fly_num} = traces_ind_fly;
end

bckg_spds = [0 10 20 40 80 160 320 640];
num_spds = length(bckg_spds);
mat_bckg_spds = repmat([0 10 20 40 80 160 0 20 40 80 160 320 0 40 80 160 320 640], num_flies,1);

mean_resp_by_bckg_spd_traces = cell(num_spds,1);
means_across_flies_resp_bckg_spd_reg_traces = zeros(25, num_spds);
sems_across_flies_resp_bckg_spd_reg_traces = zeros(25, num_spds);
for spd_idx = 1:num_spds
    curr_spd_resp_idxs = mat_bckg_spds == bckg_spds(spd_idx);
    curr_spd_resps = [];
    for fly = 1:num_flies
        curr_spd_resps = [curr_spd_resps, bckg_mean_resp_traces{fly}(:,curr_spd_resp_idxs(fly,:))];
    end
    mean_resp_by_bckg_spd_traces(spd_idx) = {curr_spd_resps};

    means_across_flies_resp_bckg_spd_reg_traces(:,spd_idx) = mean(curr_spd_resps,2);
    sems_across_flies_resp_bckg_spd_reg_traces(:,spd_idx) = std(curr_spd_resps,0,2)/sqrt(num_flies);
end


% TIME TRACE PLOTTING for BG spds
lr_bg_color = [0.7, 0.15, 0.51];
rl_bg_color = [0.11, 0.75, 0.55];
% all plots are of the first 1 second of only background motion
x_ax_time = curr_x_ax(lower_bound_idx:upper_bound_idx);
figure('Position', [1, 1, 600, 100]);
t = tiledlayout(1,num_spds,'TileSpacing','Tight');
for spd = 1:num_spds
    nexttile; 
    hold on;
    %plot bckg start 
    xline(0, "--k");
    %     xline(1000, "--k");
    ylim(time_trace_y_bounds);
    mean_spd_resp_prog_trace = means_across_flies_resp_bckg_spd_prog_traces(:,spd);
    mean_spd_resp_reg_trace = means_across_flies_resp_bckg_spd_reg_traces(:,spd);
    sem_spd_resp_prog_trace = sems_across_flies_resp_bckg_spd_prog_traces(:,spd);
    sem_spd_resp_reg_trace = sems_across_flies_resp_bckg_spd_reg_traces(:,spd);                   
    PlotXvsY(x_ax_time, mean_spd_resp_prog_trace, 'error', sem_spd_resp_prog_trace, 'color', lr_bg_color);
    PlotXvsY(x_ax_time, mean_spd_resp_reg_trace, 'error', sem_spd_resp_reg_trace, 'color', rl_bg_color);
    % axis
    %     line([0,1000],[2,2], 'Color', [0 0 0])
    %     line([2000,2000],[0,1], 'Color', [0 0 0])
    if spd == 1
        line([-500,0],[0,0], 'Color', [0 0 0])
        line([-500,-500],[0,2], 'Color', [0 0 0])
    end
    
    xlim([-500,1000])
    axis off
    hold off
end
%% FIG 5H: tuning bckg spds only
load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\sems_across_flies_resp_bckg_spd_reg");
load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\sems_across_flies_resp_bckg_spd_prog");
load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\means_resp_bckg_spd_reg");
load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\means_resp_bckg_spd_prog");
figure();
hold on
bckg_spds = [10 20 40 80 160 320 640];
errorbar(bckg_spds, mean(means_resp_bckg_spd_prog(:,2:end),1), sems_across_flies_resp_bckg_spd_prog(2:end), 'o-' , 'Color', lr_bg_color, 'CapSize',0, 'MarkerFaceColor', lr_bg_color);
errorbar(bckg_spds, mean(means_resp_bckg_spd_reg(:,2:end),1), sems_across_flies_resp_bckg_spd_reg(2:end), 'o-' , 'Color', rl_bg_color, 'CapSize',0, 'MarkerFaceColor', rl_bg_color);
errorbar([2.5], mean(means_resp_bckg_spd_prog(:,1),1), sems_across_flies_resp_bckg_spd_prog(1), 'o-' , 'Color', lr_bg_color, 'CapSize',0, 'MarkerFaceColor', lr_bg_color);
errorbar([2.5], mean(means_resp_bckg_spd_reg(:,1),1), sems_across_flies_resp_bckg_spd_reg(1), 'o-' , 'Color', rl_bg_color, 'CapSize',0, 'MarkerFaceColor', rl_bg_color);
xlabel('v_{bckg}');
ylabel('Mean Response')

legend({'v_{Background} Progressive', 'v_{Background} Regressive'})
xticks([2.5 10,20,40,80,160,320,640])
yticks(tuning_curve_yticks)
xticklabels({'0' '10' '20' '40' '80' '160' '320' '640'});
ylim(tuning_curve_y_bounds)
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
hold off

% t-test for bckg only 
p_vals_bg = zeros(1, num_spds);
for bckg_vel_idx = 1:num_spds
    [h,p,ci,stats] = ttest2(means_resp_bckg_spd_prog(:,bckg_vel_idx), means_resp_bckg_spd_reg(:,bckg_vel_idx));
    p_vals_bg(bckg_vel_idx) = p;
end

rej_null = find(p_vals_bg < 0.05);
p_vals_bg


%% FIG 5i: draw timeseries plots for rel. motion stims 

% your cell type
cellType = 'LC15';
% your indicator
sensor = 'GC6f';
% your name
surgeon = 'Elizabeth';
% your stimulus
stim = 'bar40_80_160';  
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
data_struct = RunAnalysis(args{:});
bar40_idxs = 1:6;
bar80_idxs = 7:12;
bar160_idxs = 13:18;
%get epoch info from paramfile
num_flies = size(data_struct.analysis{1, 1}.indFly,2);
param_file = readcell('C:\Users\Lab User\Documents\GitHub\psycho5\paramfiles\Elizabeth\AllByAll\bar40_80_160.txt');
epoch_names = param_file(3,3:end);
epoch_durations = param_file(7, 3:end); epoch_durations = cell2mat(epoch_durations); epoch_durations = epoch_durations./60*1000+1000;
num_epochs  = size(epoch_names, 2);


%get values for y axis
mean_resp_mat = data_struct.analysis{1, 1}.respMatPlot;
sem_resp_mat = data_struct.analysis{1, 1}.respMatSemPlot;

%make x-axis timings appropriate for bar speed
x_ax_timings = data_struct.analysis{1, 1}.timeX ;
x_ax_timings_truncated = zeros(size(x_ax_timings,1), num_epochs);
for epoch = 1:num_epochs
    timing_in_epoch = x_ax_timings;
    idx_not_in_epoch = x_ax_timings > epoch_durations(epoch);
    timing_in_epoch(idx_not_in_epoch) = NaN;
    x_ax_timings_truncated(:,epoch) = timing_in_epoch;
end

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
t = tiledlayout(3,6,'TileSpacing','Tight');
for bar40_idx = bar40_idxs
    nexttile
    hold on
    curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat(:,bar40_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat(:,bar40_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    %plot bckg start vs bar start
    %xline(0, "--", Color = [.7 .7 .7]);
    xline(1000, "--k");
    ylim(time_trace_y_bounds);
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [0.4940 0.1840 0.5560]);
    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar40_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar40_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar40_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    %plot bckg start vs bar start
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [0.4660 0.6740 0.1880]);
    %legend()
    xlim([time_trace_lower_t curr_x_ax(end)]);
    line(trace_t_scale_bar,[0 0], 'Color', [0 0 0])
    line([-2000,-2000],trace_y_scale_bar, 'Color', [0 0 0])
    axis off
    hold off
end

for bar80_idx = bar80_idxs
    nexttile
    hold on
    curr_epoch_idxs = x_ax_timings_truncated(:,bar80_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat(:,bar80_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat(:,bar80_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    %plot bckg start vs bar start
    %xline(0, "--", Color = [.7 .7 .7]);
    xline(1000, "--k");
    ylim(time_trace_y_bounds);
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [0.4940 0.1840 0.5560]);

    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar80_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar80_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar80_idx);
    curr_sem = curr_sem(curr_epoch_idxs);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    line(trace_t_scale_bar,[0 0], 'Color', [0 0 0])
    line([-2000,-2000],trace_y_scale_bar, 'Color', [0 0 0])
    axis off
    %plot bckg start vs bar start
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [0.4660 0.6740 0.1880]);
    hold off
end

for bar160_idx = bar160_idxs
    nexttile
    curr_epoch_idxs = x_ax_timings_truncated(:,bar160_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat(:,bar160_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat(:,bar160_idx);
    curr_sem = curr_sem(curr_epoch_idxs);

    %plot bckg start vs bar start
    %xline(0, "--", Color = [.7 .7 .7]);
    xline(1000, "--k");
    ylim(time_trace_y_bounds);
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [0.4940 0.1840 0.5560]);


    curr_epoch_idxs = x_ax_timings_truncated_anti(:,bar160_idx);
    curr_epoch_idxs = ~isnan(curr_epoch_idxs);

    curr_x_ax = x_ax_timings_anti;
    curr_x_ax = curr_x_ax(curr_epoch_idxs);

    curr_y_ax = mean_resp_mat_anti(:,bar160_idx);
    curr_y_ax = curr_y_ax(curr_epoch_idxs);

    curr_sem = sem_resp_mat_anti(:,bar160_idx);
    curr_sem = curr_sem(curr_epoch_idxs);
    xlim([time_trace_lower_t curr_x_ax(end)]);
    line(trace_t_scale_bar,[0 0], 'Color', [0 0 0])
    line([-2000,-2000],trace_y_scale_bar, 'Color', [0 0 0])
    axis off
    %plot bckg start vs bar start
    PlotXvsY(curr_x_ax, curr_y_ax, 'error', curr_sem, 'color', [0.4660 0.6740 0.1880]);
    hold off
end

%% FIG 5J: plot tuning curve wrt background speed, bar spd 80 only

load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\prog_relmot_means.mat");

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

%errorbar(bar40_bckgvels(2:end), (bar40_fwhms_mean(2:end)), (bar40_fwhms_sem(2:end)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[1 1 1]);
errorbar(bar80_bckgvels(2:end), (bar80_fwhms_mean(2:end)), (bar80_fwhms_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
%errorbar(bar160_bckgvels(2:end), (bar160_fwhms_ratios_mean(2:end)), (bar160_fwhms_sem(2:end)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[1 1 1]);

%errorbar(bar40_bckgvels(1), (bar40_fwhms_mean(1)), (bar40_fwhms_sem(1)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[1 1 1]);
errorbar(bar80_bckgvels(1), (bar80_fwhms_mean(1)), (bar80_fwhms_sem(1)),'o-' , 'CapSize',0, 'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
%errorbar(bar160_bckgvels(1), (bar160_fwhms_mean(1)), (bar160_fwhms_sem(1)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[1 1 1]);

all_flies_avg_resp_fwhm_opp = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\opposing_relmot_means.mat");
all_flies_avg_resp_fwhm_opp = all_flies_avg_resp_fwhm_opp.all_flies_avg_resp_fwhm;
all_flies_avg_resp_fwhm_opp = all_flies_avg_resp_fwhm_opp(1:num_flies_anti,:);

bar40_fwhms_opp = all_flies_avg_resp_fwhm_opp(:, bar40_idxs);
bar40_fwhms_mean_opp = mean(bar40_fwhms_opp, 1);
bar40_fwhms_sem_opp = std(bar40_fwhms_opp, 1)/sqrt(num_flies_anti);

bar80_fwhms_opp = all_flies_avg_resp_fwhm_opp(:, bar80_idxs);
bar80_fwhms_mean_opp = mean(bar80_fwhms_opp, 1);
bar80_fwhms_sem_opp = std(bar80_fwhms_opp, 1)/sqrt(num_flies_anti);

bar160_fwhms_opp = all_flies_avg_resp_fwhm_opp(:, bar160_idxs);
bar160_fwhms_mean_opp = mean(bar160_fwhms_opp, 1);
bar160_fwhms_sem_opp = std(bar160_fwhms_opp, 1)/sqrt(num_flies_anti);

%errorbar(bar40_bckgvels(2:end), (bar40_fwhms_mean_opp(2:end)), (bar40_fwhms_sem_opp(2:end)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[1 1 1]);
errorbar(bar80_bckgvels(2:end), (bar80_fwhms_mean_opp(2:end)), (bar80_fwhms_sem_opp(2:end)),'o-' , 'CapSize',0, 'Color', [0.4660 0.6740 0.1880], 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
%errorbar(bar160_bckgvels(2:end), (bar160_fwhms_ratios_mean(2:end)), (bar160_fwhms_sem_opp(2:end)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[1 1 1]);

%errorbar(bar40_bckgvels(1), (bar40_fwhms_mean_opp(1)), (bar40_fwhms_sem_opp(1)), 'o-' , 'CapSize',0,'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[1 1 1]);
errorbar(bar80_bckgvels(1), (bar80_fwhms_mean_opp(1)), (bar80_fwhms_sem_opp(1)),'o-' , 'CapSize',0, 'Color', [0.4660 0.6740 0.1880], 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
%errorbar(bar160_bckgvels(1), (bar160_fwhms_mean_opp(1)), (bar160_fwhms_sem_opp(1)),'o-' , 'CapSize',0, 'Color', [0.3010 0.7450 0.9330], 'MarkerFaceColor',[1 1 1]);




xlabel('v_{bckg}');
ylabel('Mean Response')
%legend({'v_{bar} = 40 °/sec', 'v_{bar} = 80 °/sec', 'v_{bar} = 160 °/sec'})
line([25,25],[0,1], 'Color', [0 0 0])
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
ylim(tuning_curve_y_bounds)
xticks([10,20,40,80,160,320,640])
yticks(tuning_curve_yticks)

% t-test for bar spd = 80
p_vals_bar_80 = zeros(1, size(bar80_bckgvels,2));
for bckg_vel_idx = 1:size(bar80_bckgvels,2)
    [h,p,ci,stats] = ttest2(bar80_fwhms(:,bckg_vel_idx), bar80_fwhms_opp(:,bckg_vel_idx));
    p_vals_bar_80(bckg_vel_idx) = p;
end
rej_null = find(p_vals_bar_80 < 0.05);
errorbar(bar80_bckgvels(rej_null), (bar80_fwhms_mean_opp(rej_null)), (bar80_fwhms_sem_opp(rej_null)),'o' , 'CapSize',0, 'Color', [0.4660 0.6740 0.1880], 'MarkerFaceColor',[0.4660 0.6740 0.1880], 'LineStyle', 'none');
p_vals_bar_80
hold off
