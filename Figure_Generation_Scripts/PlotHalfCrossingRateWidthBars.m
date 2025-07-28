controlGenotype = 3;
genotypeCounterVec = [5,1:3,40,41,12:15,17,16,18,42,7:11,19:39];

% Number of rows and columns in the tiled layout
numTilesRows = 4;
numTilesCols = ceil(length(genotypeCounterVec)/numTilesRows);

figure
TL = tiledlayout(numTilesRows,numTilesCols);

for genotypeCounter = 1:(size(genotypeCounterVec,2))
    % The line below is needed to allow user to choose order of genotypes
    genotype = genotypeCounterVec(genotypeCounter);

    crossRates_gal4_sh = zeros(length(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate),1);

    for flyCounter = 1:length(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate)
        AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:) = ...
            AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:) + 0.0001*rand(1,4);
        crossRates_gal4_sh(flyCounter) = ...
            interp1(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:), 1:0.5:2.5, 0.5);
    end
    crossRates_gal4_sh(isnan(crossRates_gal4_sh)) = 1;
    halfCrossWidth_gal4_sh = mean(crossRates_gal4_sh);

    crossRates_gal4_plus = zeros(length(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),1);

    for flyCounter = 1:length(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate)
        AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:) = ...
            AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:) + 0.0001*rand(1,4);
        crossRates_gal4_plus(flyCounter) = ...
            interp1(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:), 1:0.5:2.5, 0.5);
    end
    crossRates_gal4_plus(isnan(crossRates_gal4_plus)) = 1;
    halfCrossWidth_gal4_plus = mean(crossRates_gal4_plus);

    crossRates_empty_sh = zeros(length(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),1);

    for flyCounter = 1:length(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate)
        AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:) = ...
            AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:) + 0.0001*rand(1,4);
        crossRates_empty_sh(flyCounter) = ...
            interp1(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate(flyCounter,:), 1:0.5:2.5, 0.5);
    end
    crossRates_empty_sh(isnan(crossRates_empty_sh)) = 1;
    halfCrossWidth_empty_sh = mean(crossRates_empty_sh);

    synth_mean = ...
            mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate) + ...  % Gal4/+
            mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate) - ... % empty / sh
            mean(AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventRate); % ID1

    if sum(synth_mean == 0) > 1
            synth_mean = ...
                synth_mean + 0.0001*rand(1,4);
    end
    crossRates_synth = ...
        interp1(synth_mean, 1:0.5:2.5, 0.5);
    crossRates_synth(isnan(crossRates_synth)) = 1;
    halfCrossWidth_synth = mean(crossRates_synth);

    ax = nexttile;

    if any((5*(genotype-1)+1:5*genotype)==(p_values_most_conservative_comps_scaled<0.05).*p_values_most_conservative_comps_vec_sorted_ind,'all')
        ctrl_mean = halfCrossWidth_synth;
        ctrl_sem = std(crossRates_gal4_sh)/sqrt(length(crossRates_gal4_sh)-1);
    elseif column_of_largest_p(5*genotype) == 1 %synth
        ctrl_mean = halfCrossWidth_synth;
        ctrl_sem = std(crossRates_gal4_sh)/sqrt(length(crossRates_gal4_sh)-1);
    elseif column_of_largest_p(5*genotype) == 2 %emp>sh
        ctrl_mean = mean(crossRates_empty_sh);
        ctrl_sem = std(crossRates_empty_sh)/sqrt(length(crossRates_empty_sh)-1);
    elseif column_of_largest_p(5*genotype) == 3 %gal4>+
        ctrl_mean = mean(crossRates_gal4_plus);
        ctrl_sem = std(crossRates_gal4_plus)/sqrt(length(crossRates_gal4_plus)-1);
    end

    gal4_sh_mean = mean(crossRates_gal4_sh);
    gal4_sh_sem = std(crossRates_gal4_sh)/sqrt(length(crossRates_gal4_sh)-1);

    hold on
    data = gal4_sh_mean;
    errlow = gal4_sh_sem;
    errhigh = errlow;
    % bar(1,data,'FaceColor',[1 0 0])                
    % errorbar(1,data,errlow,errhigh,'color','k','CapSize',0);  
    barh(2,data,'FaceColor',[1 0 0])                
    errorbar(data,2,0,0,errlow,errhigh,'color','k','CapSize',0);  
    % er.Color = [0 0 0];                 
    % er.LineStyle = 'none';   
    data = ctrl_mean;
    errlow = ctrl_sem;
    errhigh = errlow;
    % bar(2,data,'FaceColor',[.75 .75 .75]) 
    % errorbar(2,data,errlow,errhigh,'Color','k','CapSize',0);  
    barh(1,data,'FaceColor',[.75 .75 .75]) 
    errorbar(data,1,0,0,errlow,errhigh,'Color','k','CapSize',0);  
    % er.Color = [0 0 0];                            
    % er.LineStyle = 'none';
    title(erase(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)),'split '));
    hold off
    xlim([0,2.5])
    ylim([0.5,2.5])
    yticks([]);
    if genotypeCounter+numTilesCols >= length(genotypeCounterVec)
        xticks([0,2.5])
    else
        xticks([]);
    end

end

TL.XLabel.String = 'Width of Half Crossing Probability (mm)';