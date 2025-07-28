controlGenotype = 3;
genotypeCounterVec = [5,1:3,40,41,12:15,17,16,18,42,7:11,19:39];

% Number of rows and columns in the tiled layout
numTilesRows = 4;
numTilesCols = ceil(length(genotypeCounterVec)/numTilesRows);

figure
TL = tiledlayout(numTilesRows,numTilesCols);

% Go through each individual genotype now and grab the data and perform
% statistical tests
for genotypeCounter = 1:(size(genotypeCounterVec,2))
    % The line below is needed to allow user to choose order of genotypes
    genotype = genotypeCounterVec(genotypeCounter);

    synth_mean = ...
            mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate,'all') + ...  % Gal4/+
            mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate,'all') - ... % empty / sh
            mean(AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventRate,'all'); % ID1
    synth_sem = sqrt(...
            sum((AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
            (AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
            (AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;
    
    gal4_plus_mean = ...
            mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate,'all');
    gal4_plus_sem = sqrt(...
            sum((AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;
    
    empty_sh_mean = ...
            mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate,'all');
    empty_sh_sem = sqrt(...
            sum((AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;
    
    gal4_sh_mean = ...
            mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate,'all');
    gal4_sh_sem = sqrt(...
            sum((AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;

    ax = nexttile;

    if any((5*(genotype-1)+1:5*genotype)==(p_values_most_conservative_comps_scaled<0.05).*p_values_most_conservative_comps_vec_sorted_ind,'all')
        ctrl_mean = synth_mean;
        ctrl_sem = synth_sem;
    elseif column_of_largest_p(5*genotype) == 1 %synth
        ctrl_mean = synth_mean;
        ctrl_sem = synth_sem;
    elseif column_of_largest_p(5*genotype) == 2 %emp>sh
        ctrl_mean = empty_sh_mean;
        ctrl_sem = empty_sh_sem;
    elseif column_of_largest_p(5*genotype) == 3 %gal4>+
        ctrl_mean = gal4_plus_mean;
        ctrl_sem = gal4_plus_sem;
    end

    hold on
    data = gal4_sh_mean;
    errlow = gal4_sh_sem;
    errhigh = errlow;
    bar(1,data,'FaceColor',[1 0 0])                
    errorbar(1,data,errlow,errhigh,'color','k');  
    % er.Color = [0 0 0];                 
    % er.LineStyle = 'none';   
    data = ctrl_mean;
    errlow = ctrl_sem;
    errhigh = errlow;
    bar(2,data,'FaceColor',[.75 .75 .75]) 
    errorbar(2,data,errlow,errhigh,'Color','k');  
    % er.Color = [0 0 0];                            
    % er.LineStyle = 'none';
    title(erase(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)),'split '));
    hold off
    ylim([0,.6])
    xlim([0.6,2.4])
    xticks([]);
    if mod(genotypeCounter-1,numTilesCols) == 0
        yticks([0,0.6])
        ylabel('Crossing Curve Average')
    else
        yticks([]);
    end

end

TL.TileSpacing = 'compact';
TL.Padding = 'compact';

figure
TL2 = tiledlayout(numTilesRows,numTilesCols);
% gapSizeVec = [1:0.5:2.5]';
gapSizeVec = [1,1,1,1]';
% Go through each individual genotype now and grab the data and perform
% statistical tests
for genotypeCounter = 1:(size(genotypeCounterVec,2))
    hit = 0;
    % The line below is needed to allow user to choose order of genotypes
    genotype = genotypeCounterVec(genotypeCounter);

    synth_mean = ...
            (mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate)*gapSizeVec + ...  % Gal4/+
            mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate)*gapSizeVec - ... % empty / sh
            mean(AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventRate)*gapSizeVec)/sum(gapSizeVec); % ID1
    synth_sem = sqrt(...
            sum((AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
            (AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
            (AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;
    
    gal4_plus_mean = ...
            mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate)*gapSizeVec/sum(gapSizeVec);
    gal4_plus_sem = sqrt(...
            sum((AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;
    
    empty_sh_mean = ...
            mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate)*gapSizeVec/sum(gapSizeVec);
    empty_sh_sem = sqrt(...
            sum((AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;
    
    gal4_sh_mean = ...
            mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate)*gapSizeVec/sum(gapSizeVec);
    gal4_sh_sem = sqrt(...
            sum((AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd).^2))/4;

    ax = nexttile;

    if any((5*(genotype-1)+1:5*genotype)==(p_values_most_conservative_comps_scaled<0.05).*p_values_most_conservative_comps_vec_sorted_ind,'all')
        ctrl_mean = synth_mean;
        ctrl_sem = synth_sem;
        hit = 1;
    elseif column_of_largest_p(5*genotype) == 1 %synth
        ctrl_mean = synth_mean;
        ctrl_sem = synth_sem;
    elseif column_of_largest_p(5*genotype) == 2 %emp>sh
        ctrl_mean = empty_sh_mean;
        ctrl_sem = empty_sh_sem;
    elseif column_of_largest_p(5*genotype) == 3 %gal4>+
        ctrl_mean = gal4_plus_mean;
        ctrl_sem = gal4_plus_sem;
    end

    hold on
    data = gal4_sh_mean;
    errlow = gal4_sh_sem;
    errhigh = errlow;
    bar(1,data,'FaceColor',[1 0 0])                
    errorbar(1,data,errlow,errhigh,'color','k','CapSize',0);  
    % er.Color = [0 0 0];                 
    % er.LineStyle = 'none';   
    data = ctrl_mean;
    errlow = ctrl_sem;
    errhigh = errlow;
    bar(2,data,'FaceColor',[.8 .8 .8]) 
    errorbar(2,data,errlow,errhigh,'Color','k','CapSize',0);  
    % er.Color = [0 0 0];                            
    % er.LineStyle = 'none';
    xlabel(erase(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)),'split '));
    if hit
        % plot([1,2],[0.5,0.5],'k');
        plot([1,2],[0.75,0.75],'k');
    end
    hold off
    ylim([0,1])
    xlim([0.5,2.5])
    xticks([]);
    if mod(genotypeCounter-1,numTilesCols) == 0
        % yticks([0,0.6])
        % ylabel('Crossing Curve CoM')
        yticks([0,1])
        ylabel('Average Crossing Probability')
    else
        yticks([]);
        ax.YAxis.Visible = 'off';
    end

    set(gca, 'color', 'none');
    % set(gca,'visible','off')
    ax = gca;
    set(gca,'xtick',[])
    

end

TL2.TileSpacing = 'compact';
TL2.Padding = 'compact';