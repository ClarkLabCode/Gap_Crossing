%% FIGS6E: compare features of quiescent
compFeaturesRsq("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\quiescent_flies_avg_resp_fwhm.mat")
%% FIGS6F: compare features of active
compFeaturesRsq("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\active_flies_avg_resp_fwhm.mat")
%%
function [ratio_bar, ratio_bg, bg_bar] = compFeaturesRsq(meanFileDir)

    %compare r_sq predictivity of stimulus features (ratios of bar speed+background speed, bar speed, background speed) used to predict neural response 
    all_flies_avg_resp_fwhm = load(meanFileDir);
    fname = fieldnames(all_flies_avg_resp_fwhm);
    all_flies_avg_resp_fwhm = all_flies_avg_resp_fwhm.(fname{1});
    num_flies = size(all_flies_avg_resp_fwhm,1);
    
    %normalize within fly
    all_flies_avg_resp_fwhm = (all_flies_avg_resp_fwhm./max(all_flies_avg_resp_fwhm,[],2 ));
    all_flies_avg_resp_fwhm = all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    %first calculate observed rsq vals by stim feature type: this is the empirical
    % RATIOS
    v_ratio_x_ax = log([0.25 0.5 1 2 4]);
    y_ratios_bckg_bar = repmat([v_ratio_x_ax v_ratio_x_ax v_ratio_x_ax], num_flies,1);
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    
    %BIN BY NEURAL RESPONSE
    
    %get idxs of each bin
    resp_idxs_by_bin = zeros([size(all_flies_avg_resp_fwhm), num_bins]); %dims: [flies, conditions, bins]
    for bin_idx = 1:num_bins
        if bin_idx == 1
            resp_idxs_by_bin(:,:,bin_idx) = all_flies_avg_resp_fwhm <= resp_bins(bin_idx);
        else
            resp_idxs_by_bin(:,:,bin_idx) = logical((all_flies_avg_resp_fwhm <= resp_bins(bin_idx)) .* (all_flies_avg_resp_fwhm > resp_bins(bin_idx-1)));
        end
    
    end
    
    %get neural resp by bin
    neural_resp_by_bin = cell(num_bins,1);
    for bin_idx = 1:num_bins
        neural_resp_by_bin(bin_idx) = {all_flies_avg_resp_fwhm(logical(resp_idxs_by_bin(:,:,bin_idx)))};
    
    end
    
    %get speed ratios corresponding to nerual resp bin
    spd_ratios_by_bin = cell(num_bins,1);
    for bin_idx = 1:num_bins
        spd_ratios_by_bin(bin_idx) = {y_ratios_bckg_bar(logical(resp_idxs_by_bin(:,:,bin_idx)))};
    
    end
    
    
    %nonparametric r^2 calculation
    y_bar = mean(mean(y_ratios_bckg_bar));
    sq_res_arr = [];
    sq_tot_arr = [];
    
    %calculate square residuals and square total for each bin of neural response
    for bin_idx = 1:num_bins
        y_hat_bin = mean(spd_ratios_by_bin{bin_idx}); %MLE
    
        %iter thru each point in neural resp bin
        for yi=spd_ratios_by_bin{bin_idx}'
            sq_res = (yi-y_hat_bin)^2; %for SSR in numerator
            sq_tot = (yi-y_bar)^2; %for SSE in denom
            sq_res_arr = [sq_res_arr sq_res];
            sq_tot_arr = [sq_tot_arr sq_tot];
            
        end
    end
    r_sq_ratio_obs = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
    % plot decoder
    
    y_ratio_bckg_bar_means = zeros(1,num_bins);
    y_ratio_bckg_bar_sems = zeros(1,num_bins);
    
    for bin_idx = 1:num_bins
        y_ratio_bckg_bar_means(bin_idx) = mean(spd_ratios_by_bin{bin_idx});
        y_ratio_bckg_bar_sems(bin_idx) = std(spd_ratios_by_bin{bin_idx})/sqrt(num_flies);
    end
    
    figure()
    hold on
    errorbar(resp_bins-0.05, y_ratio_bckg_bar_means, y_ratio_bckg_bar_sems, 'o-' , 'CapSize',0,'Color', [0 0 0], 'MarkerFaceColor',[0 0 0], 'LineStyle', 'none');
    title(strcat('r^{2}= ',num2str(r_sq_ratio_obs)))
    xlim([0 1.1]);
    xlabel('Neural response bin');
    xticks(resp_bins)
    %xticklabels({'0 - .1','.1 - .2','.2 - .3','.4 - .5','.6 - .7','.7 - .8','.8 - .9','.9 - 1'})
    yticks(v_ratio_x_ax)
    ylabel('v_{background}/v_{bar}');
    yticklabels({'0.25','0.5','1','2', '4'})
    hold off
    
    % BAR SPEEDS
    bar_spds = [40 80 160];
    y_bar_spds = repmat([40 40 40 40 40 80 80 80 80 80 160 160 160 160 160], num_flies,1);
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    %BIN BY NEURAL RESPONSE
    
    %get idxs of each bin
    resp_idxs_by_bin = zeros([size(all_flies_avg_resp_fwhm), num_bins]); %dims: [flies, conditions, bins]
    for bin_idx = 1:num_bins
        if bin_idx == 1
            resp_idxs_by_bin(:,:,bin_idx) = all_flies_avg_resp_fwhm <= resp_bins(bin_idx);
        else
            resp_idxs_by_bin(:,:,bin_idx) = logical((all_flies_avg_resp_fwhm <= resp_bins(bin_idx)) .* (all_flies_avg_resp_fwhm > resp_bins(bin_idx-1)));
        end
    
    end
    
    %get neural resp by bin
    neural_resp_by_bin = cell(num_bins,1);
    for bin_idx = 1:num_bins
        neural_resp_by_bin(bin_idx) = {all_flies_avg_resp_fwhm(logical(resp_idxs_by_bin(:,:,bin_idx)))};
    
    end
    
    %get bar spds corresponding to nerual resp bin
    bar_spd_by_bin = cell(num_bins,1);
    for bin_idx = 1:num_bins
        bar_spd_by_bin(bin_idx) = {y_bar_spds(logical(resp_idxs_by_bin(:,:,bin_idx)))};
    end
    
    
    
    %nonparametric r^2 calculation
    y_bar = mean(mean(y_bar_spds));
    sq_res_arr = [];
    sq_tot_arr = [];
    
    %calculate square residuals and square total for each bin of neural response
    for bin_idx = 1:num_bins
        y_hat_bin = mean(bar_spd_by_bin{bin_idx}); %MLE
    
        %iter thru each point in neural resp bin
        for yi=bar_spd_by_bin{bin_idx}'
            sq_res = (yi-y_hat_bin)^2; %for SSR in numerator
            sq_tot = (yi-y_bar)^2; %for SSE in denom
            sq_res_arr = [sq_res_arr sq_res];
            sq_tot_arr = [sq_tot_arr sq_tot];
            
        end
    end
    r_sq_bar_spd_obs = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
    % plot decoder
    
    y_bar_spd_means = zeros(1,num_bins);
    y_bar_spd_sems = zeros(1,num_bins);
    
    for bin_idx = 1:num_bins
        y_bar_spd_means(bin_idx) = mean(bar_spd_by_bin{bin_idx});
        y_bar_spd_sems(bin_idx) = std(bar_spd_by_bin{bin_idx})/sqrt(num_flies);
    end
    
    figure()
    hold on
    errorbar(resp_bins-0.05, y_bar_spd_means, y_bar_spd_sems, 'o-' , 'CapSize',0,'Color', [0 0 0], 'MarkerFaceColor',[0 0 0], 'LineStyle', 'none');
    set(gca, 'YScale', 'log')
    title(strcat('r^{2}= ',num2str(r_sq_bar_spd_obs)))
    xlim([0 1.1]);
    xlabel('Neural response bin');
    xticks(resp_bins)
    %xticklabels({'0 - .1','.1 - .2','.2 - .3','.4 - .5','.6 - .7','.7 - .8','.8 - .9','.9 - 1'})
    yticks(bar_spds)
    ylim([20 320])
    ylabel('v_{bar}');
    yticklabels({'40','80','160'})
    hold off
    
    % BCKG SPD
    bckg_spds = [10 20 40 80 160 320 640];
    y_bckg_spds = repmat([10 20 40 80 160 20 40 80 160 320 40 80 160 320 640], num_flies,1);
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    %BIN BY NEURAL RESPONSE
    
    %get idxs of each bin
    resp_idxs_by_bin = zeros([size(all_flies_avg_resp_fwhm), num_bins]); %dims: [flies, conditions, bins]
    for bin_idx = 1:num_bins
        if bin_idx == 1
            resp_idxs_by_bin(:,:,bin_idx) = all_flies_avg_resp_fwhm <= resp_bins(bin_idx);
        else
            resp_idxs_by_bin(:,:,bin_idx) = logical((all_flies_avg_resp_fwhm <= resp_bins(bin_idx)) .* (all_flies_avg_resp_fwhm > resp_bins(bin_idx-1)));
        end
    
    end
    
    %get neural resp by bin
    neural_resp_by_bin = cell(num_bins,1);
    for bin_idx = 1:num_bins
        neural_resp_by_bin(bin_idx) = {all_flies_avg_resp_fwhm(logical(resp_idxs_by_bin(:,:,bin_idx)))};
    
    end
    
    %get bar spds corresponding to nerual resp bin
    bckg_spd_by_bin = cell(num_bins,1);
    for bin_idx = 1:num_bins
        bckg_spd_by_bin(bin_idx) = {y_bckg_spds(logical(resp_idxs_by_bin(:,:,bin_idx)))};
    
    end
    
    
    %nonparametric r^2 calculation
    y_bar = mean(mean(y_bckg_spds));
    sq_res_arr = [];
    sq_tot_arr = [];
    
    
    
    %calculate square residuals and square total for each bin of neural response
    for bin_idx = 1:num_bins
        y_hat_bin = mean(bckg_spd_by_bin{bin_idx}); %MLE
    
        %iter thru each point in neural resp bin
        for yi=bckg_spd_by_bin{bin_idx}'
            sq_res = (yi-y_hat_bin)^2; %for SSR in numerator
            sq_tot = (yi-y_bar)^2; %for SSE in denom
            sq_res_arr = [sq_res_arr sq_res];
            sq_tot_arr = [sq_tot_arr sq_tot];
            
        end
    end
    r_sq_bckg_spd_obs = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
    %plot decoder
    
    y_bckg_spd_means = zeros(1,num_bins);
    y_bckg_spd_sems = zeros(1,num_bins);
    
    for bin_idx = 1:num_bins
        y_bckg_spd_means(bin_idx) = mean(bckg_spd_by_bin{bin_idx});
        y_bckg_spd_sems(bin_idx) = std(bckg_spd_by_bin{bin_idx})/sqrt(length(bckg_spd_by_bin{bin_idx}));
    end
    
    figure()
    hold on
    errorbar(resp_bins-0.05, y_bckg_spd_means, y_bckg_spd_sems, 'o-' , 'CapSize',0,'Color', [0 0 0], 'MarkerFaceColor',[0 0 0], 'LineStyle', 'none');
    set(gca, 'YScale', 'log')
    title(strcat('r^{2}= ',num2str(r_sq_bckg_spd_obs)))
    xlim([0 1.1]);
    xlabel('Neural response bin');
    xticks(resp_bins)
    %xticklabels({'0 - .1','.1 - .2','.2 - .3','.4 - .5','.6 - .7','.7 - .8','.8 - .9','.9 - 1'})
    yticks(bckg_spds)
    ylim([0 1280])
    ylabel('v_{background}');
    
    hold off
    %bar spd
    bar_spds = [40 80 160];
    y_bar_spds = repmat([40 40 40 40 40 80 80 80 80 80 160 160 160 160 160], num_flies,1);
    
    bckg_spds = [10 20 40 80 160 320 640];
    y_bckg_spds = repmat([10 20 40 80 160 20 40 80 160 320 40 80 160 320 640], num_flies,1);
    
    v_ratio_x_ax = log([0.25 0.5 1 2 4]);     % calc wrt. speed ratios, note in log space so we exclude cases where condition = 0 
    y_ratios_bckg_bar = repmat([v_ratio_x_ax v_ratio_x_ax v_ratio_x_ax], num_flies,1);
    
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    n_iter = 10000;
    r_boostrapped_bar_spds = zeros(1, n_iter);
    r_boostrapped_bckg_spds = zeros(1, n_iter);
    r_boostrapped_ratios = zeros(1, n_iter);
    for iter = 1:n_iter
        
        bootstrap_idxs = randi(num_flies,[num_flies,1]);
        curr_iter_all_flies_avg_resp_fwhm = all_flies_avg_resp_fwhm(bootstrap_idxs,:);
    
        
        %BIN BY NEURAL RESPONSE    
        %get idxs of each bin
        resp_idxs_by_bin = zeros([size(curr_iter_all_flies_avg_resp_fwhm), num_bins]); %dims: [flies, conditions, bins]
        for bin_idx = 1:num_bins
            if bin_idx == 1
                resp_idxs_by_bin(:,:,bin_idx) = curr_iter_all_flies_avg_resp_fwhm <= resp_bins(bin_idx);
            else
                resp_idxs_by_bin(:,:,bin_idx) = logical((curr_iter_all_flies_avg_resp_fwhm <= resp_bins(bin_idx)) .* (curr_iter_all_flies_avg_resp_fwhm > resp_bins(bin_idx-1)));
            end
        
        end
        
        %get neural resp by bin
        neural_resp_by_bin = cell(num_bins,1);
        for bin_idx = 1:num_bins
            neural_resp_by_bin(bin_idx) = {curr_iter_all_flies_avg_resp_fwhm(logical(resp_idxs_by_bin(:,:,bin_idx)))};
        
        end
        
        %get bar spds corresponding to nerual resp bin
        bar_spd_by_bin = cell(num_bins,1);
        for bin_idx = 1:num_bins
            bar_spd_by_bin(bin_idx) = {y_bar_spds(logical(resp_idxs_by_bin(:,:,bin_idx)))};
        end
        
        
    
        %nonparametric r^2 calculation
        y_bar = mean(mean(y_bar_spds));
        sq_res_arr = [];
        sq_tot_arr = [];
        
        
        
        %calculate square residuals and square total for each bin of neural response
        for bin_idx = 1:num_bins
            y_hat_bin = mean(bar_spd_by_bin{bin_idx}); %MLE
        
            %iter thru each point in neural resp bin
            for yi=bar_spd_by_bin{bin_idx}'
                sq_res = (yi-y_hat_bin)^2; %for SSR in numerator
                sq_tot = (yi-y_bar)^2; %for SSE in denom
                sq_res_arr = [sq_res_arr sq_res];
                sq_tot_arr = [sq_tot_arr sq_tot];
                
            end
        end
        r_sq_bar_spd = 1-sum(sq_res_arr)/sum(sq_tot_arr);
        r_boostrapped_bar_spds(iter) = r_sq_bar_spd ;
    
        %get bg spds corresponding to neural resp bin
        bckg_spd_by_bin = cell(num_bins,1);
        for bin_idx = 1:num_bins
            bckg_spd_by_bin(bin_idx) = {y_bckg_spds(logical(resp_idxs_by_bin(:,:,bin_idx)))};
        
        end
    
        %nonparametric r^2 calculation
        y_bar = mean(mean(y_bckg_spds));
        sq_res_arr = [];
        sq_tot_arr = [];
        
        
        
        %calculate square residuals and square total for each bin of neural response
        for bin_idx = 1:num_bins
            y_hat_bin = mean(bckg_spd_by_bin{bin_idx}); %MLE
        
            %iter thru each point in neural resp bin
            for yi=bckg_spd_by_bin{bin_idx}'
                sq_res = (yi-y_hat_bin)^2; %for SSR in numerator
                sq_tot = (yi-y_bar)^2; %for SSE in denom
                sq_res_arr = [sq_res_arr sq_res];
                sq_tot_arr = [sq_tot_arr sq_tot];
                
            end
        end
        r_sq_bckg_spd = 1-sum(sq_res_arr)/sum(sq_tot_arr);
        r_boostrapped_bckg_spds(iter) = r_sq_bckg_spd ;
    
    
         %get speed ratios corresponding to nerual resp bin
        spd_ratios_by_bin = cell(num_bins,1);
        for bin_idx = 1:num_bins
            spd_ratios_by_bin(bin_idx) = {y_ratios_bckg_bar(logical(resp_idxs_by_bin(:,:,bin_idx)))};
        
        end
    
        
        %nonparametric r^2 calculation
        y_bar = mean(mean(y_ratios_bckg_bar));
        sq_res_arr = [];
        sq_tot_arr = [];
        
        %calculate square residuals and square total for each bin of neural response
        for bin_idx = 1:num_bins
            y_hat_bin = mean(spd_ratios_by_bin{bin_idx}); %MLE
        
            %iter thru each point in neural resp bin
            for yi=spd_ratios_by_bin{bin_idx}'
                sq_res = (yi-y_hat_bin)^2; %for SSR in numerator
                sq_tot = (yi-y_bar)^2; %for SSE in denom
                sq_res_arr = [sq_res_arr sq_res];
                sq_tot_arr = [sq_tot_arr sq_tot];
                
            end
        end
        r_sq_ratio = 1-sum(sq_res_arr)/sum(sq_tot_arr);
        r_boostrapped_ratios(iter) = r_sq_ratio ;
    end
    
    %get CIs and plot histss
    
    figure()
    hold on
    histogram(r_boostrapped_bar_spds);
    title("CI R^{2} of Bar Spd")
    xline(r_sq_bar_spd_obs, Color='red' );
    bars_ci = prctile(r_boostrapped_bar_spds, [2.5 97.5])
    xline(prctile(r_boostrapped_bar_spds, [2.5]), Color='black' );
    xline(prctile(r_boostrapped_bar_spds, [97.5]), Color='black' );
    xlabel("");
    hold off
    
    figure()
    hold on
    histogram(r_boostrapped_bckg_spds);
    title("CI R^{2} of BG Spd")
    xline(r_sq_bckg_spd_obs, Color='red' );
    bckg_ci = prctile(r_boostrapped_bckg_spds, [2.5 97.5])
    xline(prctile(r_boostrapped_bckg_spds, [2.5]), Color='black' );
    xline(prctile(r_boostrapped_bckg_spds, [97.5]), Color='black' );
    xlabel("");
    hold off
    
    figure()
    hold on
    histogram(r_boostrapped_ratios);
    title("CI R^{2} Ratios")
    xline(r_sq_ratio_obs, Color='red' );
    ratios_ci = prctile(r_boostrapped_ratios, [2.5 97.5])
    xline(prctile(r_boostrapped_ratios, [2.5]), Color='black' );
    xline(prctile(r_boostrapped_ratios, [97.5]), Color='black' );
    xlabel("");
    hold off
    %get p-vals
    
    
    ratio_bar_obs = r_sq_ratio_obs - r_sq_bar_spd_obs;
    ratio_bar = sum(abs(r_boostrapped_ratios - r_boostrapped_bar_spds) >   abs(2*ratio_bar_obs))/n_iter
    
    ratio_bg_obs = r_sq_ratio_obs - r_sq_bckg_spd_obs;
    ratio_bg = sum(abs(r_boostrapped_ratios - r_boostrapped_bckg_spds) >   abs(2*ratio_bg_obs))/n_iter
    
    bg_bar_obs = r_sq_bckg_spd_obs - r_sq_bar_spd_obs;
    bg_bar = sum(abs(r_boostrapped_bckg_spds - r_boostrapped_bar_spds) >   abs(2*bg_bar_obs))/n_iter
    
    figure()
    hold on
    rs = [r_sq_bar_spd_obs r_sq_bckg_spd_obs r_sq_ratio_obs];
    errorbar([1 2 3], rs, abs(rs-[bars_ci(1) bckg_ci(1) ratios_ci(1)]), abs(rs-[bars_ci(2) bckg_ci(2) ratios_ci(2)]), 'o' , 'CapSize',0,'Color', [0 0 0], 'MarkerFaceColor',[0 0 0])
    yticks([0, 0.25, 0.5]);
    xticks([1 2 3])
    xticklabels({'Bar speed', 'BG Speed', 'Ratio'})
    title('Feature R^{2}')
    xlim([0.5 3.5])
    hold off
end


