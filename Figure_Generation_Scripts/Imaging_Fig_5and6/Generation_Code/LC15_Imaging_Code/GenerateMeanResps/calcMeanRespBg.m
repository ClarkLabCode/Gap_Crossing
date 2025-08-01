%% 
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
%dataPath = 'Y:\2p_microscope_data\w_+;LC15AD_GC6f;LC15DBD_+\bar40_80_160\2023\10_18\12_42_55';
%roiExtractionFile = 'WatershedRegionRestrictedRoiExtraction'; % use watershed + manual circling around
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
%% 
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
lower_bound_msec = 0; %only background displayed from 0 to 1 sec
upper_bound_msec = 1000;
    
%loop through flies, then loops through timeseries of response to each
%epoch. Takes only the first 1000 msec of stim presentation because that is
%the time in which only the bckg is moving (bar has not started yet)
bckg_mean_resp = zeros(num_flies, num_epochs);
for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p8_averagedRois.snipMat{epoch, 1};
    end
    

    means_ind_fly = zeros(1, num_epochs);

    
    for bar40_idx = bar40_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
    

    
        %integrate
        area = trapz(resp_window)*acq_rate; 
        delta_time_sec = (curr_x_ax(upper_bound_idx) - curr_x_ax(lower_bound_idx))/1000;
        mean_value_resp = area/delta_time_sec;
    
        means_ind_fly(bar40_idx) = mean_value_resp;
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
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
    

    
        %integrate
        area = trapz(resp_window)*acq_rate; 
        delta_time_sec = (curr_x_ax(upper_bound_idx) - curr_x_ax(lower_bound_idx))/1000;
        mean_value_resp = area/delta_time_sec;
    
        means_ind_fly(bar80_idx) = mean_value_resp;
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
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
    

    
        %integrate
        area = trapz(resp_window)*acq_rate; 
        delta_time_sec = (curr_x_ax(upper_bound_idx) - curr_x_ax(lower_bound_idx))/1000;
        mean_value_resp = area/delta_time_sec;
    
        means_ind_fly(bar160_idx) = mean_value_resp;
    end
    


    bckg_mean_resp(fly_num,:) = means_ind_fly;
end

%% prog backgrounds
bckg_spds = [0 10 20 40 80 160 320 640];
num_spds = length(bckg_spds);
%mat_bckg_spds = repmat([0 10 20 40 80 160 0 20 40 80 160 320 0 40 80 160 320 640], num_flies,1);
mat_bckg_spds = [0 10 20 40 80 160 0 20 40 80 160 320 0 40 80 160 320 640];


means_resp_bckg_spd_prog = zeros(num_flies, num_spds);

%since bckg spds get repeated per experiment we should combine them within
%flies
for fly_num = 1:num_flies
    curr_fly_mean_resps = bckg_mean_resp(fly_num,:);
    for spd_idx = 1:num_spds
        curr_spd_resp_idxs = mat_bckg_spds == bckg_spds(spd_idx);
        curr_spd_resps = curr_fly_mean_resps(curr_spd_resp_idxs);
        means_resp_bckg_spd_prog(fly_num,spd_idx) = mean(curr_spd_resps);
    

    end


end

sems_across_flies_resp_bckg_spd_prog = std(means_resp_bckg_spd_prog,1)/sqrt(num_flies);

%% reg

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

% get path to the data
dataPath = GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
%dataPath = 'Y:\2p_microscope_data\w_+;LC15AD_GC6f;LC15DBD_+\bar40_80_160\2023\10_18\12_42_55';
%roiExtractionFile = 'WatershedRegionRestrictedRoiExtraction'; % use watershed + manual circling around
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
lower_bound_msec = 0; %only background displayed from 0 to 1 sec
upper_bound_msec = 1000;
    
%loop through flies, then loops through timeseries of response to each
%epoch. Takes only the first 1000 msec of stim presentation because that is
%the time in which only the bckg is moving (bar has not started yet)
bckg_mean_resp = zeros(num_flies, num_epochs);
for fly_num = 1:num_flies
    ind_mean_resp_mat = zeros(size(mean_resp_mat));
    for epoch = 1:num_epochs
        ind_mean_resp_mat(:,epoch) = a1.analysis{1, 1}.indFly{1, fly_num}.p8_averagedRois.snipMat{epoch, 1};
    end
    

    means_ind_fly = zeros(1, num_epochs);

    
    for bar40_idx = bar40_idxs
        curr_epoch_idxs = x_ax_timings_truncated(:,bar40_idx);
        curr_epoch_idxs = ~isnan(curr_epoch_idxs);
    
        curr_x_ax = x_ax_timings;
        curr_x_ax = curr_x_ax(curr_epoch_idxs);
    
        curr_y_ax = ind_mean_resp_mat(:,bar40_idx);
        curr_y_ax = curr_y_ax(curr_epoch_idxs);
    
        %convert input bounds (in msec) to sampling rate that directly
        %correspond to resp mat idxs
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
    

    
        %integrate
        area = trapz(resp_window)*acq_rate; 
        delta_time_sec = (curr_x_ax(upper_bound_idx) - curr_x_ax(lower_bound_idx))/1000;
        mean_value_resp = area/delta_time_sec;
    
        means_ind_fly(bar40_idx) = mean_value_resp;
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
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
    

    
        %integrate
        area = trapz(resp_window)*acq_rate; 
        delta_time_sec = (curr_x_ax(upper_bound_idx) - curr_x_ax(lower_bound_idx))/1000;
        mean_value_resp = area/delta_time_sec;
    
        means_ind_fly(bar80_idx) = mean_value_resp;
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
        lower_bound_idx = sum(curr_x_ax < lower_bound_msec);
        upper_bound_idx = sum(curr_x_ax < upper_bound_msec);
    
        resp_window = curr_y_ax(lower_bound_idx:upper_bound_idx);
    

    
        %integrate
        area = trapz(resp_window)*acq_rate; 
        delta_time_sec = (curr_x_ax(upper_bound_idx) - curr_x_ax(lower_bound_idx))/1000;
        mean_value_resp = area/delta_time_sec;
    
        means_ind_fly(bar160_idx) = mean_value_resp;
    end
    


    bckg_mean_resp(fly_num,:) = means_ind_fly;
end

%% reg backgrounds
bckg_spds = [0 10 20 40 80 160 320 640];
num_spds = length(bckg_spds);
%mat_bckg_spds = repmat([0 10 20 40 80 160 0 20 40 80 160 320 0 40 80 160 320 640], num_flies,1);
mat_bckg_spds = [0 10 20 40 80 160 0 20 40 80 160 320 0 40 80 160 320 640];

means_resp_bckg_spd_reg = zeros(num_flies, num_spds);

%since bckg spds get repeated per experiment we should combine them within
%flies
for fly_num = 1:num_flies
    curr_fly_mean_resps = bckg_mean_resp(fly_num,:);
    for spd_idx = 1:num_spds
        curr_spd_resp_idxs = mat_bckg_spds == bckg_spds(spd_idx);
        curr_spd_resps = curr_fly_mean_resps(curr_spd_resp_idxs);
        means_resp_bckg_spd_reg(fly_num,spd_idx) = mean(curr_spd_resps);
    

    end


end

sems_across_flies_resp_bckg_spd_reg = std(means_resp_bckg_spd_reg,1)/sqrt(num_flies);

%%

figure();
hold on
errorbar([2.5 10 20 40 80 160 320 640], mean(means_resp_bckg_spd_prog,1), sems_across_flies_resp_bckg_spd_prog, 'o-' , 'Color', [0 0.4470 0.7410], 'CapSize',0, 'MarkerFaceColor', [0 0.4470 0.7410]);
errorbar([2.5 10 20 40 80 160 320 640], mean(means_resp_bckg_spd_reg,1), sems_across_flies_resp_bckg_spd_reg, 'o-' , 'Color', [0.8 0.2470 0.15], 'CapSize',0, 'MarkerFaceColor', [0.8 0.2470 0.15]);
xlabel('v_{bckg}');
ylabel('Mean Response')

legend({'v_{Background} Progressive', 'v_{Background} Regressive'})
xticks([2.5 10,20,40,80,160,320,640])
xticklabels({'0' '10' '20' '40' '80' '160' '320' '640'});
ylim([-0.5 4.5])
set(gca, 'XScale', 'log')
hold off

p_vals = zeros(1, num_spds);
for bckg_vel_idx = 1:num_spds
    [h,p,ci,stats] = ttest2(means_resp_bckg_spd_prog(:,bckg_vel_idx), means_resp_bckg_spd_reg(:,bckg_vel_idx));
    p_vals(bckg_vel_idx) = p;
end

rej_null = find(p_vals < 0.05);
p_vals
