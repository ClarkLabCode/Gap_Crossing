    %%  r_sq predictivity of the ratio of bar:background speed used to predict neural response, compared across ftbftb, ftb btf, and btf btf regimes
    %% ftb ftb
    load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\prog_relmot_means.mat");
    
    
    %normalize within fly
    all_flies_avg_resp_fwhm = (all_flies_avg_resp_fwhm./max(all_flies_avg_resp_fwhm,[],2 ));
    all_flies_avg_resp_fwhm = all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies = size(all_flies_avg_resp_fwhm,1);
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
    
    %get speed ratios corresponding to neural resp bin
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
    r_sq_ratio_obs_ftbftb = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
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
    title(strcat('r^{2}= ',num2str(r_sq_ratio_obs_ftbftb)))
    xlim([0 1.1]);
    xlabel('Neural response bin');
    xticks(resp_bins)
    %xticklabels({'0 - .1','.1 - .2','.2 - .3','.4 - .5','.6 - .7','.7 - .8','.8 - .9','.9 - 1'})
    yticks(v_ratio_x_ax)
    ylabel('v_{background}/v_{bar}');
    yticklabels({'0.25','0.5','1','2', '4'})
    hold off
    
    
    n_iter = 10000;
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
        
        %get speed ratios corresponding to neural resp bin
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
    
    %plot hist
    figure()
    hold on
    histogram(r_boostrapped_ratios);
    title("Bootstrap for R^{2} of Ratios")
    xline(r_sq_ratio_obs_ftbftb, Color='red' );
    ratios_ci_ftbftb = prctile(r_boostrapped_ratios, [2.5 97.5])
    xline(prctile(r_boostrapped_ratios, [2.5]), Color='black' );
    xline(prctile(r_boostrapped_ratios, [97.5]), Color='black' );
    xlabel("");
    hold off
    
    %% btf ftb
    load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\opposing_relmot_means.mat");
    
    
    %normalize within fly
    all_flies_avg_resp_fwhm = (all_flies_avg_resp_fwhm./max(all_flies_avg_resp_fwhm,[],2 ));
    all_flies_avg_resp_fwhm = all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies = size(all_flies_avg_resp_fwhm,1);
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
    
    %get speed ratios corresponding to neural resp bin
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
    r_sq_ratio_obs_ftbbtf = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
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
    title(strcat('r^{2}= ',num2str(r_sq_ratio_obs_ftbbtf)))
    xlim([0 1.1]);
    xlabel('Neural response bin');
    xticks(resp_bins)
    %xticklabels({'0 - .1','.1 - .2','.2 - .3','.4 - .5','.6 - .7','.7 - .8','.8 - .9','.9 - 1'})
    yticks(v_ratio_x_ax)
    ylabel('v_{background}/v_{bar}');
    yticklabels({'0.25','0.5','1','2', '4'})
    hold off
    
    n_iter = 10000;
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
        
        %get speed ratios corresponding to neural resp bin
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
    
    %plot hist
    figure()
    hold on
    histogram(r_boostrapped_ratios);
    title("Bootstrap for R^{2} of Ratios")
    xline(r_sq_ratio_obs_ftbbtf, Color='red' );
    ratios_ci_ftbbtf = prctile(r_boostrapped_ratios, [2.5 97.5])
    xline(prctile(r_boostrapped_ratios, [2.5]), Color='black' );
    xline(prctile(r_boostrapped_ratios, [97.5]), Color='black' );
    xlabel("");
    hold off
    
    
    %% btf btf
    load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\flipped_relmot_means.mat");
    num_flies = size(all_flies_avg_resp_fwhm,1);
    
    %normalize within fly
    all_flies_avg_resp_fwhm = (all_flies_avg_resp_fwhm./max(all_flies_avg_resp_fwhm,[],2 ));
    all_flies_avg_resp_fwhm = all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    
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
    
    %get speed ratios corresponding to neural resp bin
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
    r_sq_ratio_obs_btfbtf = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
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
    title(strcat('r^{2}= ',num2str(r_sq_ratio_obs_btfbtf)))
    xlim([0 1.1]);
    xlabel('Neural response bin');
    xticks(resp_bins)
    %xticklabels({'0 - .1','.1 - .2','.2 - .3','.4 - .5','.6 - .7','.7 - .8','.8 - .9','.9 - 1'})
    yticks(v_ratio_x_ax)
    ylabel('v_{background}/v_{bar}');
    yticklabels({'0.25','0.5','1','2', '4'})
    hold off
    
    n_iter = 10000;
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
        
        %get speed ratios corresponding to neural resp bin
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
    
    %plot hist
    figure()
    hold on
    histogram(r_boostrapped_ratios);
    title("Bootstrap for R^{2} of Ratios")
    xline(r_sq_ratio_obs_btfbtf, Color='red' );
    ratios_ci_btfbtf = prctile(r_boostrapped_ratios, [2.5 97.5])
    xline(prctile(r_boostrapped_ratios, [2.5]), Color='black' );
    xline(prctile(r_boostrapped_ratios, [97.5]), Color='black' );
    xlabel("");
    hold off
    
    
    %%
    figure()
    hold on
    rs = [r_sq_ratio_obs_ftbftb r_sq_ratio_obs_ftbbtf r_sq_ratio_obs_btfbtf];
    errorbar([1 2 3], rs, abs(rs-[ratios_ci_ftbftb(1) ratios_ci_ftbbtf(1) ratios_ci_btfbtf(1)]), abs(rs-[ratios_ci_ftbftb(2) ratios_ci_ftbbtf(2) ratios_ci_btfbtf(2)]), 'o' , 'CapSize',0,'Color', [0 0 0], 'MarkerFaceColor',[0 0 0])
    yticks([0, 0.25, 0.5]);
    xticklabels({'FTB-FTB', 'FTB-BTF', 'BTF-BTF'})
    xticks([1 2 3])
    xlim([0.5 3.5])
    hold off
    
    
    
    
    %% ftb ftb minus btf ftb
    
    %ftbftb
    ftb_ftb_all_flies_avg_resp_fwhm = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\prog_relmot_means.mat");
    ftb_ftb_all_flies_avg_resp_fwhm = ftb_ftb_all_flies_avg_resp_fwhm.all_flies_avg_resp_fwhm;
    %normalize within fly
    ftb_ftb_all_flies_avg_resp_fwhm = (ftb_ftb_all_flies_avg_resp_fwhm./max(ftb_ftb_all_flies_avg_resp_fwhm,[],2 ));
    ftb_ftb_all_flies_avg_resp_fwhm = ftb_ftb_all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies_regime1 = size(ftb_ftb_all_flies_avg_resp_fwhm,1);
    
    %ftbbtf
    ftb_btf_all_flies_avg_resp_fwhm = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\opposing_relmot_means.mat");
    ftb_btf_all_flies_avg_resp_fwhm = ftb_btf_all_flies_avg_resp_fwhm.all_flies_avg_resp_fwhm;
    %normalize within fly
    ftb_btf_all_flies_avg_resp_fwhm = (ftb_btf_all_flies_avg_resp_fwhm./max(ftb_btf_all_flies_avg_resp_fwhm,[],2 ));
    ftb_btf_all_flies_avg_resp_fwhm = ftb_btf_all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies_regime2 = size(ftb_btf_all_flies_avg_resp_fwhm,1);
    
    %pooling for null distribution
    pooled_all_flies_avg_resp_fwhm = [ftb_ftb_all_flies_avg_resp_fwhm; ftb_btf_all_flies_avg_resp_fwhm];
    num_flies_pooled = num_flies_regime1+num_flies_regime2;
    
    v_ratio_x_ax = log([0.25 0.5 1 2 4]);
    y_ratios_bckg_bar = repmat([v_ratio_x_ax v_ratio_x_ax v_ratio_x_ax], num_flies_pooled,1);
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    
    n_iter = 10000;
    rsq_diffs_bootstrapped = zeros(1, n_iter);
    
    for iter = 1:n_iter
        bootstrap_idxs = randi(num_flies_pooled,[num_flies_pooled,1]);
        curr_iter_all_flies_avg_resp_fwhm_pooled = pooled_all_flies_avg_resp_fwhm(bootstrap_idxs,:);
    
        %CALCULATIONS FOR REGIME 1
        curr_iter_all_flies_avg_resp_fwhm = curr_iter_all_flies_avg_resp_fwhm_pooled(1:num_flies_regime1,:);
    
        %bin by neural resp
        
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
        
        %get speed ratios corresponding to neural resp bin
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
        r_sq_regime1 = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
        
        %CALCULATIONS FOR REGIME 2
        curr_iter_all_flies_avg_resp_fwhm = curr_iter_all_flies_avg_resp_fwhm_pooled((num_flies_regime1+1):num_flies_pooled,:);
    
        %bin by neural resp
        
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
        
        %get speed ratios corresponding to neural resp bin
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
        r_sq_regime2 = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
        r_sq_diff_regimes = r_sq_regime1-r_sq_regime2;
        rsq_diffs_bootstrapped(iter) = r_sq_diff_regimes;
    end
    
    
    %two tailed
    empirical_diff = r_sq_ratio_obs_ftbftb - r_sq_ratio_obs_ftbbtf;
    p_val_ftbftb_minus_ftbbtf = 1-sum(abs(empirical_diff) >= abs(rsq_diffs_bootstrapped))/n_iter
    
    %% ftb ftb minus btf btf
    
    %ftbftb
    ftb_ftb_all_flies_avg_resp_fwhm = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\prog_relmot_means.mat");
    ftb_ftb_all_flies_avg_resp_fwhm = ftb_ftb_all_flies_avg_resp_fwhm.all_flies_avg_resp_fwhm;
    %normalize within fly
    ftb_ftb_all_flies_avg_resp_fwhm = (ftb_ftb_all_flies_avg_resp_fwhm./max(ftb_ftb_all_flies_avg_resp_fwhm,[],2 ));
    ftb_ftb_all_flies_avg_resp_fwhm = ftb_ftb_all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies_regime1 = size(ftb_ftb_all_flies_avg_resp_fwhm,1);
    
    %btfbtf
    btf_btf_all_flies_avg_resp_fwhm = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\flipped_relmot_means.mat");
    btf_btf_all_flies_avg_resp_fwhm = btf_btf_all_flies_avg_resp_fwhm.all_flies_avg_resp_fwhm;
    %normalize within fly
    btf_btf_all_flies_avg_resp_fwhm = (btf_btf_all_flies_avg_resp_fwhm./max(btf_btf_all_flies_avg_resp_fwhm,[],2 ));
    btf_btf_all_flies_avg_resp_fwhm = btf_btf_all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies_regime2 = size(btf_btf_all_flies_avg_resp_fwhm,1);
    
    %pooling for null distribution
    pooled_all_flies_avg_resp_fwhm = [ftb_ftb_all_flies_avg_resp_fwhm; btf_btf_all_flies_avg_resp_fwhm];
    num_flies_pooled = num_flies_regime1+num_flies_regime2;
    
    v_ratio_x_ax = log([0.25 0.5 1 2 4]);
    y_ratios_bckg_bar = repmat([v_ratio_x_ax v_ratio_x_ax v_ratio_x_ax], num_flies_pooled,1);
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    
    n_iter = 10000;
    rsq_diffs_bootstrapped = zeros(1, n_iter);
    
    for iter = 1:n_iter
        bootstrap_idxs = randi(num_flies_pooled,[num_flies_pooled,1]);
        curr_iter_all_flies_avg_resp_fwhm_pooled = pooled_all_flies_avg_resp_fwhm(bootstrap_idxs,:);
    
        %CALCULATIONS FOR REGIME 1
        curr_iter_all_flies_avg_resp_fwhm = curr_iter_all_flies_avg_resp_fwhm_pooled(1:num_flies_regime1,:);
    
        %bin by neural resp
        
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
        
        %get speed ratios corresponding to neural resp bin
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
        r_sq_regime1 = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
        
        %CALCULATIONS FOR REGIME 2
        curr_iter_all_flies_avg_resp_fwhm = curr_iter_all_flies_avg_resp_fwhm_pooled((num_flies_regime1+1):num_flies_pooled,:);
    
        %bin by neural resp
        
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
        
        %get speed ratios corresponding to neural resp bin
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
        r_sq_regime2 = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
        r_sq_diff_regimes = r_sq_regime1-r_sq_regime2;
        rsq_diffs_bootstrapped(iter) = r_sq_diff_regimes;
    end
    
    
    %two tailed
    empirical_diff = r_sq_ratio_obs_ftbftb - r_sq_ratio_obs_btfbtf;
    p_val_ftbftb_minus_btfbtf = 1-sum(abs(empirical_diff) >= abs(rsq_diffs_bootstrapped))/n_iter
    
    
    
        
    %% ftb btf minus btf btf
    
    %ftbbtf
    ftb_btf_all_flies_avg_resp_fwhm = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\opposing_relmot_means.mat");
    ftb_btf_all_flies_avg_resp_fwhm = ftb_btf_all_flies_avg_resp_fwhm.all_flies_avg_resp_fwhm;
    %normalize within fly
    ftb_btf_all_flies_avg_resp_fwhm = (ftb_btf_all_flies_avg_resp_fwhm./max(ftb_btf_all_flies_avg_resp_fwhm,[],2 ));
    ftb_btf_all_flies_avg_resp_fwhm = ftb_btf_all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies_regime1 = size(ftb_btf_all_flies_avg_resp_fwhm,1);
    
    %btfbtf
    btf_btf_all_flies_avg_resp_fwhm = load("C:\Users\Lab User\Documents\GitHub\psycho5\analysis\analysisFiles\Elizabeth\mean_resps_data_structs\flipped_relmot_means.mat");
    btf_btf_all_flies_avg_resp_fwhm = btf_btf_all_flies_avg_resp_fwhm.all_flies_avg_resp_fwhm;
    %normalize within fly
    btf_btf_all_flies_avg_resp_fwhm = (btf_btf_all_flies_avg_resp_fwhm./max(btf_btf_all_flies_avg_resp_fwhm,[],2 ));
    btf_btf_all_flies_avg_resp_fwhm = btf_btf_all_flies_avg_resp_fwhm(:, [2:6 8:12 14:18]);
    num_flies_regime2 = size(btf_btf_all_flies_avg_resp_fwhm,1);
    
    %pooling for null distribution
    pooled_all_flies_avg_resp_fwhm = [ftb_btf_all_flies_avg_resp_fwhm; btf_btf_all_flies_avg_resp_fwhm];
    num_flies_pooled = num_flies_regime1+num_flies_regime2;
    
    v_ratio_x_ax = log([0.25 0.5 1 2 4]);
    y_ratios_bckg_bar = repmat([v_ratio_x_ax v_ratio_x_ax v_ratio_x_ax], num_flies_pooled,1);
    resp_bins = [0.1:.1:1]; %10 bins
    num_bins = length(resp_bins);
    
    
    n_iter = 10000;
    rsq_diffs_bootstrapped = zeros(1, n_iter);
    
    for iter = 1:n_iter
        bootstrap_idxs = randi(num_flies_pooled,[num_flies_pooled,1]);
        curr_iter_all_flies_avg_resp_fwhm_pooled = pooled_all_flies_avg_resp_fwhm(bootstrap_idxs,:);
    
        %CALCULATIONS FOR REGIME 1
        curr_iter_all_flies_avg_resp_fwhm = curr_iter_all_flies_avg_resp_fwhm_pooled(1:num_flies_regime1,:);
    
        %bin by neural resp
        
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
        
        %get speed ratios corresponding to neural resp bin
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
        r_sq_regime1 = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
        
        %CALCULATIONS FOR REGIME 2
        curr_iter_all_flies_avg_resp_fwhm = curr_iter_all_flies_avg_resp_fwhm_pooled((num_flies_regime1+1):num_flies_pooled,:);
    
        %bin by neural resp
        
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
        
        %get speed ratios corresponding to neural resp bin
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
        r_sq_regime2 = 1-sum(sq_res_arr)/sum(sq_tot_arr);
    
        r_sq_diff_regimes = r_sq_regime1-r_sq_regime2;
        rsq_diffs_bootstrapped(iter) = r_sq_diff_regimes;
    end
    
    
    %two tailed
    empirical_diff = r_sq_ratio_obs_ftbbtf - r_sq_ratio_obs_btfbtf;
    p_val_ftbbtf_minus_btfbtf = 1-sum(abs(empirical_diff) >= abs(rsq_diffs_bootstrapped))/n_iter

