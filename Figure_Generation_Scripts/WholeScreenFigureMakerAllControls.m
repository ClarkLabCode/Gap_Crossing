% WHOLESCREENFIGUREMAKER Creates the figure of all genotypes' crossing curves
%
%  This script generates the main screen figure for the paper. It is a tiled
%  figure in which each tile corresponds to a different neuron being silenced
%  along with its control(s). The user is asked to select which control(s)
%  to show in the figure (the options are the synthetic control, both
%  parental controls, or the most statistically conservative control).
%
%  In order to successfully run this script, the cell arrays AllCrossStats,
%  AllCrossStatsNames, and AllCrossStatsGenotypes must all be loaded into 
%  the workspace.
%
%  This function also generates all the statistical info relevant for the
%  screen. To see the Bonferroni-Holm corrected statistical tests, look at
%  the following two output variables:
%     p_values_most_conservative_comps_scaled       <-- Bonferroni-Holm corrected p values
%     p_values_most_conservative_comps_sorted_names <-- Corresponding neuron and comparison

% Check that the user has loaded in AllCrossStats
if ~exist('AllCrossStats','var')
    error('Please load in AllCrossStats, AllCrossStatsNames, and AllCrossStatsGenotypes before running this script.')
end

% Initialize vectors that will hold p values for point-by-point comparisons
p_of_points_to_synth_vec_up = ones(4*(max(size(AllCrossStats))),1);
p_of_points_to_empty_over_shts_vec_up = ones(4*(max(size(AllCrossStats))),1);
p_of_points_to_gal4_over_plus_vec_up = ones(4*(max(size(AllCrossStats))),1);

% Initialize vectors that will hold p values for curve-to-curve comparisons
p_of_curve_diff_from_synth_vec_up = ones((max(size(AllCrossStats))),1);
p_of_curve_diff_from_empty_over_shts_vec_up = ones((max(size(AllCrossStats))),1);
p_of_curve_diff_from_gal4_over_plus_vec_up = ones((max(size(AllCrossStats))),1);

% Allows user to choose which genotypes to plot and in what order
% genotypeCounterVec = [1:3,5,12:18,42,7:11,19:41];
genotypeCounterVec = [5,1:3,40,41,12:15,17,16,18,42,7:11,19:39];

% Number of rows and columns in the tiled layout
numTilesRows = 8;
numTilesCols = 5;

% Go through each individual genotype now and grab the data and perform
% statistical tests
for genotypeCounter = 1:(size(genotypeCounterVec,2))
    % The line below is needed to allow user to choose order of genotypes
    genotype = genotypeCounterVec(genotypeCounter);
    % If there's no data for a given selected genotype over shts, skip it
    if isempty(AllCrossStatsGenotypes{genotype,2})
        continue
    % If there's no data for a given selected genotype over IsoD1, skip it
    elseif isempty(AllCrossStatsGenotypes{genotype,3})
        continue
    % Use split empty > shts as control for split Gal4s
    elseif contains(AllCrossStatsNames{genotype,2}, 'split')
        controlGenotype = 3;
    % Use split empty DBD > shts as control for single Gal4s
    else
        controlGenotype = 4;
    end
    
    % Construct the synthetic controls (this is a model where the effect of
    % each construct [UAS and Gal4] have an individual additive effect)
    % The mean of the synthetic control algebraically simplifies to <Gal4/+> + <UAS/+> - <+> 
    synth_mean = ...
        mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate) + ...
        mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate) - ...
        mean(AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventRate);
    % Ensure that the crossing probability is non-negative
    synth_mean = synth_mean.*(synth_mean>=0);
    % The error of the synthetic control is computed via error propogation
    synth_sem = sqrt(...
        (AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
        (AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
        (AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventStd).^2);

    % Initialize a vector that holds the p values of comparing the silenced
    % experiment to the two parental controls point-by-point
    p_of_points_to_empty_over_shts_up = zeros(4,1);
    p_of_points_to_gal4_over_plus_up = zeros(4,1);
    
    % Initialize a vector that holds the z scores of comparing the silenced
    % experiment to the synthetic control point-by-point
    z_of_points_to_synth_up = zeros(4,1);

    % Go through each gap now and perform the statistical tests
    for gapSize = 1:4
        % Wilcoxon rank sum test when comparing to sp empty (DBD) > shts
        p_of_points_to_empty_over_shts_up(gapSize) = ...
            ranksum(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate(:,gapSize),...
                    AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate(:,gapSize));
        % Wilcoxon rank sum test when comparing to Gal4 > +
        p_of_points_to_gal4_over_plus_up(gapSize) = ...
            ranksum(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate(:,gapSize),...
                    AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate(:,gapSize));
        % z test when comparing to synthetic control
        z_of_points_to_synth_up(gapSize) = ...
            (mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate(:,gapSize)) - ...
             synth_mean(gapSize))./synth_sem(gapSize);
    end
    
    % Convert the z score to a p value for the comparison to synthetic control
    p_of_points_to_synth_up = 2*(1-normcdf(abs(z_of_points_to_synth_up)));

    % Now fill in the appropriate entries in the vectors that hold all 
    % p values for all of the different genotypes
    p_of_points_to_synth_vec_up((4*(genotype-1)+1):(4*genotype)) = p_of_points_to_synth_up';
    p_of_points_to_empty_over_shts_vec_up((4*(genotype-1)+1):(4*genotype)) = p_of_points_to_empty_over_shts_up';
    p_of_points_to_gal4_over_plus_vec_up((4*(genotype-1)+1):(4*genotype)) = p_of_points_to_gal4_over_plus_up';
    
    % Now we compute the curve-to-curve difference by computing the
    % cumulative difference between the curves and propogating error of
    % each curve
    % First do this for the synthetic curve-to-curve comparison
    curve_diff_from_synth = ...
        sum(mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate)) - ...
        sum(synth_mean);
    error_of_curve_diff_from_synth = sqrt(...
        sum(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd.^2) + ...
        sum(synth_sem.^2));
    % Compute the z score then convert to p value and fill in the appropriate vector
    z_score_of_curve_diff_from_synth = curve_diff_from_synth/error_of_curve_diff_from_synth;
    p_of_curve_diff_from_synth = 2*(1-normcdf(abs(z_score_of_curve_diff_from_synth)));
    p_of_curve_diff_from_synth_vec_up(genotype) = p_of_curve_diff_from_synth; 

    % Next do this for the sp empty (DBD) > shts curve-to-curve comparison
    curve_diff_from_empty_over_shts = ...
        sum(mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate)) - ...
        sum(mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate));
    error_of_curve_diff_from_empty_over_shts = sqrt(...
        sum(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd.^2) + ...
        sum(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd.^2));
    % Compute the z score then convert to p value and fill in the appropriate vector
    z_score_of_curve_diff_from_empty_over_shts = curve_diff_from_empty_over_shts/error_of_curve_diff_from_empty_over_shts;
    p_of_curve_diff_from_empty_over_shts = 2*(1-normcdf(abs(z_score_of_curve_diff_from_empty_over_shts)));
    p_of_curve_diff_from_empty_over_shts_vec_up(genotype) = p_of_curve_diff_from_empty_over_shts;    
    
    % Next do this for the Gal4 > + curve-to-curve comparison
    curve_diff_from_gal4_over_plus = ...
        sum(mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate)) - ...
        sum(mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate));
    error_of_curve_diff_from_gal4_over_plus = sqrt(...
        sum(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd.^2) + ...
        sum(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd.^2));
    % Compute the z score then convert to p value and fill in the appropriate vector
    z_score_of_curve_diff_from_gal4_over_plus = curve_diff_from_gal4_over_plus/error_of_curve_diff_from_gal4_over_plus;
    p_of_curve_diff_from_gal4_over_plus = 2*(1-normcdf(abs(z_score_of_curve_diff_from_gal4_over_plus)));
    p_of_curve_diff_from_gal4_over_plus_vec_up(genotype) = p_of_curve_diff_from_gal4_over_plus;    
    
end

% Now we combine the p values of the two different types of comparisons
% (point-to-point and curve-to-curve) into one vector of p values 
% Here we do it for the synthetic control
p_values_all_comps_synth = ones(size([p_of_points_to_synth_vec_up; p_of_curve_diff_from_synth_vec_up]));
genoCounter = 1;
% The structure of this vector is:
% [(Genotype 1: Gaps 1-4, Whole Curve), (Genotype 2: Gaps 1-4, Whole Curve),...]
for i = 1:(5*max(size(AllCrossStats)))
    if mod(i,5) ~= 0
        p_values_all_comps_synth(i) = p_of_points_to_synth_vec_up(mod(i,5)+4*(genoCounter-1));
    else
        p_values_all_comps_synth(i) = p_of_curve_diff_from_synth_vec_up(i/5);
        genoCounter = genoCounter + 1;
    end
end

% Now we combine the p values of the two different types of comparisons
% (point-to-point and curve-to-curve) into one vector of p values 
% Here we do it for the sp empty (DBD) > shts control
p_values_all_comps_empty_over_shts = ones(size([p_of_points_to_empty_over_shts_vec_up; p_of_curve_diff_from_empty_over_shts_vec_up]));
genoCounter = 1;
% The structure of this vector is:
% [(Genotype 1: Gaps 1-4, Whole Curve), (Genotype 2: Gaps 1-4, Whole Curve),...]
for i = 1:(5*max(size(AllCrossStats)))
    if mod(i,5) ~= 0
        p_values_all_comps_empty_over_shts(i) = p_of_points_to_empty_over_shts_vec_up(mod(i,5)+4*(genoCounter-1));
    else
        p_values_all_comps_empty_over_shts(i) = p_of_curve_diff_from_empty_over_shts_vec_up(i/5);
        genoCounter = genoCounter + 1;
    end
end

% Now we combine the p values of the two different types of comparisons
% (point-to-point and curve-to-curve) into one vector of p values 
% Here we do it for the Gal4 > + control
p_values_all_comps_gal4_over_plus = ones(size([p_of_points_to_gal4_over_plus_vec_up; p_of_curve_diff_from_gal4_over_plus_vec_up]));
genoCounter = 1;
% The structure of this vector is:
% [(Genotype 1: Gaps 1-4, Whole Curve), (Genotype 2: Gaps 1-4, Whole Curve),...]
for i = 1:(5*max(size(AllCrossStats)))
    if mod(i,5) ~= 0
        p_values_all_comps_gal4_over_plus(i) = p_of_points_to_gal4_over_plus_vec_up(mod(i,5)+4*(genoCounter-1));
    else
        p_values_all_comps_gal4_over_plus(i) = p_of_curve_diff_from_gal4_over_plus_vec_up(i/5);
        genoCounter = genoCounter + 1;
    end
end

% Now we make a matrix that holds all of these comparisons to the 3 controls
% (col 1 = synth, col 2 = sp empty (DBD) > shts, col 3 = Gal4 > +)
p_values_all_comps = [p_values_all_comps_synth, p_values_all_comps_empty_over_shts, p_values_all_comps_gal4_over_plus];

% We want to choose the comparison that is the most conservative, so look
% for which of the three comparisons has the highest minimum p value
column_of_largest_p = zeros(size(p_values_all_comps,1)/5,1);
for i = 1:(size(p_values_all_comps,1)/5)
    [~,column_of_largest_p(i)] = max([min(p_values_all_comps(((5*(i-1)+1):5*i),1)),min(p_values_all_comps(((5*(i-1)+1):5*i),2)),min(p_values_all_comps(((5*(i-1)+1):5*i),3))]);
end

% For the sake of easily using this info later, we make a vector to hold
% the info of which comparison was the most conservative and we make it the
% same size as the vector that holds all p values
column_of_largest_p = repmat(column_of_largest_p,[1,5])';
column_of_largest_p = reshape(column_of_largest_p,[size(p_values_all_comps,1),1]);

% Now generate a vector of all the p values for each genotype that uses the
% most conservative comparison
p_values_most_conservative_comps_vec = zeros(size(p_values_all_comps,1),1);
for i = 1:size(p_values_all_comps,1)
    p_values_most_conservative_comps_vec(i) = p_values_all_comps(i,column_of_largest_p(i));
end

% Now we apply Holm-Bonferroni multiple comparisons correction
% To do this, we must take all our p values, sort them from smallest to
% largest, then multiply each by [the number of comparisons made minus that
% comparison's rank in p value from smallest to largest]
[p_values_most_conservative_comps_vec_sorted, p_values_most_conservative_comps_vec_sorted_ind] = sort(p_values_most_conservative_comps_vec);
p_values_most_conservative_comps_scaled = p_values_most_conservative_comps_vec_sorted.*...
                            ((length(p_values_most_conservative_comps_vec_sorted)+1-(1:length(p_values_most_conservative_comps_vec_sorted)))');

% Now that we have everything ranked and Holm-Bonferroni corrected, we
% create a cell array that holds the name of which comparison each p value
% corresponds to within the sorted and scaled list
p_values_most_conservative_comps_sorted_names = cell(size(p_values_most_conservative_comps_vec_sorted));
for i = 1:length(p_values_most_conservative_comps_vec_sorted)
    % For a given genotype, entries 1-4 correspond to Gaps 1-4
    if mod(p_values_most_conservative_comps_vec_sorted_ind(i),5) ~= 0
        p_values_most_conservative_comps_sorted_names{i} = ...
            [(AllCrossStatsNames{ceil(p_values_most_conservative_comps_vec_sorted_ind(i)/5),2}...
             (1:(strfind(AllCrossStatsNames{ceil(p_values_most_conservative_comps_vec_sorted_ind(i)/5),2}, '>')-2))),...
              ' Gap ',num2str(mod(p_values_most_conservative_comps_vec_sorted_ind(i),5))];
    % For a given genotype, entry correspond to the whole curve
    else
        p_values_most_conservative_comps_sorted_names{i} = ...
            [(AllCrossStatsNames{ceil(p_values_most_conservative_comps_vec_sorted_ind(i)/5),2}...
             (1:(strfind(AllCrossStatsNames{ceil(p_values_most_conservative_comps_vec_sorted_ind(i)/5),2}, '>')-2))),...
              ' Whole Curve'];
    end
end

% For ease of plotting this in the rest of the function below, make a
% vector that holds which comparison was made for each genotype
p_values_most_conservative_comps_sorted_control_column = column_of_largest_p(p_values_most_conservative_comps_vec_sorted_ind);

% Now that we have all the statistics finished running, go through all the
% data and plot it in the tiled layout
for genotypeCounter = 1:(size(genotypeCounterVec,2))
    % Again, this line below makes this compatible with the user's desired
    % order of genotypes to display
    genotype = genotypeCounterVec(genotypeCounter);
    % Skip genotypes if the data from silencing experiment isn't present 
    % but still put the title to maintain the desired grid shape
    if isempty(AllCrossStatsGenotypes{genotype,2})
        nexttile
        ylim([0,1]);
        xlim([0.8,2.7]);
        xticks(1:0.5:2.5);
        title(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)));
        continue
    % Skip genotypes if the data from the Gal4 > + experiment isn't present 
    % but still put the title to maintain the desired grid shape
    elseif isempty(AllCrossStatsGenotypes{genotype,3})
        nexttile
        ylim([0,1]);
        xlim([0.8,2.7]);
        xticks(1:0.5:2.5);
        title(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)));
        continue
    % Use split empty > shts as control for split Gal4s
    elseif contains(AllCrossStatsNames{genotype,2}, 'split')
        controlGenotype = 3;
    % Use split empty DBD > shts as control for single Gal4s
    else
        controlGenotype = 4;
    end
        
    % Construct the synthetic controls (this is a model where the effect of
    % each construct [UAS and Gal4] have an individual additive effect)
    % The mean of the synthetic control algebraically simplifies to <Gal4/+> + <UAS/+> - <+> 
    synth_mean = ...
        mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate) + ...
        mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate) - ...
        mean(AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventRate);
    % Ensure that the crossing probability is non-negative
    synth_mean = synth_mean.*(synth_mean>=0);
    % The error of the synthetic control is computed via error propogation
    synth_sem = sqrt(...
        (AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
        (AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
        (AllCrossStats{1,1}.AllUpVecProperCrossOverAllGapEventStd).^2);

    % Give the user the option to choose between the 3 plotting types:
    % 1) Always plot the most conservative curve (or synth if p < 0.05)
    % 2) Always plot the synthetic control curve
    % 3) Always plot the two genetic control curves
    if genotypeCounter == 1
        whichControlsToPlot = questdlg('What control curve(s) would you like plotted?', ...
	        'Choice of plotted control(s)', ...
	        'Most conservative','Synthetic','All','All');

        % Create 8x5 tiled layout to hold each curve
        figure
        % TL = tiledlayout(8,5);
        TL = tiledlayout(numTilesRows,numTilesCols);
    end

    % Move forward in the tiled layout
    ax = nexttile;

    % Parse the response
    switch whichControlsToPlot
        % If most conservative was chosen, then plot the most conservative
        % curve unless all curves are statistically different, in which
        % case we plot the synthetic curve (reasoning explained below)
        case 'Most conservative'

            % To get a global legend for the tiled layout, plot an invisible set of
            % curves that correspond to the correct color scheme of the curves
            if genotypeCounter == 1
                hold on
                    h1 = ...
                        plot(1,nan(1,1),'Color','r','Marker','none','LineWidth',5);
                    h2 = ...
                        plot(1,nan(1,1),'Color','k','Marker','none','LineWidth',5);
                    lg  = legend([h1 h2],{'Silenced','Control'},'AutoUpdate','off');
                lg.FontSize = 10;
                lg.Layout.Tile = 'North'; % <-- Legend placement with tiled layout
                lg.Orientation = 'horizontal';
                hold off
            end

            % If the silencing experiment was statistically different than both
            % parental controls and the synthetic control, plot the synthetic
            % control as the control (because the statistical significance rules  
            % out a dominant/recessive inheritance hypothesis, so parental controls 
            % are not the appropriate control to show)

            hold on
            if any((5*(genotype-1)+1:5*genotype)==(p_values_most_conservative_comps_scaled<0.05).*p_values_most_conservative_comps_vec_sorted_ind,'all')
                errorbar(1:0.5:2.5,...
                     synth_mean,...
                     synth_sem,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
                plot(1:0.5:2.5,...
                        synth_mean,'k')
            % If the silencing experiment was most conservatively compared to the
            % synthetic control, plot synthetic control as the control
            elseif column_of_largest_p(5*genotype) == 1
                errorbar(1:0.5:2.5,...
                     synth_mean,...
                     synth_sem,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
                plot(1:0.5:2.5,...
                        synth_mean,'k')
            % If the silencing experiment was most conservatively compared to the
            % sp empty (DBD) > shts control, plot sp empty (DBD) > shts as the control
            elseif column_of_largest_p(5*genotype) == 2
                errorbar(1:0.5:2.5,...
                     mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
                plot(1:0.5:2.5,...
                        mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),'k')
            % If the silencing experiment was most conservatively compared to the
            % Gal4 > + control, plot Gal4 > + as the control
            elseif column_of_largest_p(5*genotype) == 3
                errorbar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
                plot(1:0.5:2.5,...
                        mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),'k')
            end

        % If synthetic chosen, then plot the synth curve and legend
        case 'Synthetic'

            % To get a global legend for the tiled layout, plot an invisible set of
            % curves that correspond to the correct color scheme of the curves
            if genotypeCounter == 1
                hold on
                    h1 = ...
                        plot(1,nan(1,1),'Color','r','Marker','none','LineWidth',5);
                    h2 = ...
                        plot(1,nan(1,1),'Color','k','Marker','none','LineWidth',5);
                    lg  = legend([h1 h2],{'Silenced','Synthetic Control'},'AutoUpdate','off');
                lg.FontSize = 10;
                lg.Layout.Tile = 'North'; % <-- Legend placement with tiled layout
                lg.Orientation = 'horizontal';
                hold off
            end

            % Synthetic control
            hold on
            errorbar(1:0.5:2.5,...
                     synth_mean,...
                     synth_sem,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
            plot(1:0.5:2.5,...
                     synth_mean,'k')

        % If both genetic lines chosen, then plot both and legend
        case 'All'

            % To get a global legend for the tiled layout, plot an invisible set of
            % curves that correspond to the correct color scheme of the curves 
            if genotypeCounter == 1
                hold on
                h1 = ...
                    plot(1,nan(1,1),'Color','r','Marker','none','LineWidth',5);
                h2 = ...
                    plot(1,nan(1,1),'Color',[.7 .7 .7],'Marker','none','LineWidth',5);
                h3 = ...
                    plot(1,nan(1,1),'Color','k','Marker','none','LineWidth',5);
                h4 = ...
                    plot(1,nan(1,1),'k--','Marker','none','LineWidth',5);
                lg  = legend([h1 h2 h3 h4],{'[Neuron] > +;sh^{ts};shR', ...
                             '[Neuron] > +;+;+',...
                             'split empty (DBD) > +;sh^{ts};shR',...
                             'synthetic control'},'AutoUpdate','off'); 
                lg.FontSize = 10;
                        lg.Layout.Tile = 'North'; % <-- Legend placement with tiled layout
                        lg.Orientation = 'horizontal';
                hold off
            end

            % Gal4 > + control (gray)
            hold on
            errorbar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd,...
                     'Color',[.7 .7 .7],'CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
            plot(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),'Color',[.7 .7 .7])
            % sp empty (DBD) > shts control (solid black)
            errorbar(1:0.5:2.5,...
                     mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
            plot(1:0.5:2.5,...
                     mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),'k')
            % Synthetic control (dashed black)
            errorbar(1:0.5:2.5,...
                     synth_mean,...
                     synth_sem,...
                     'k','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
            plot(1:0.5:2.5,...
                     synth_mean,'k--')
    end

    % Plot the silencing experiment in red on the same plot as control(s)
    errorbar(1:0.5:2.5,...
             mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate),...
             AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd,...
             'r','CapSize',0,'LineWidth',3,'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','LineStyle','none');
    plot(1:0.5:2.5,...
             mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate),'r')
    hold off

    % Just adjust the axes and add the appropriate title now for each tile
    % genotypeCounter
    ylim([0,1]);
    xlim([0.8,2.7]);
    if genotypeCounter > (numTilesRows-1)*numTilesCols
        xticks(1:0.5:2.5);
    else
        xticks([]);
    end
    if mod(genotypeCounter,numTilesCols) == 1
        yticks([0,1]);
    else
        yticks([]);
    end
    title(erase(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)),'split '));

end

% Add a global x and y axis title since all experiments share those axes
TL.XLabel.String = 'Gap Width (mm)';
TL.YLabel.String = 'Proper Crossings / All Gap Events';
TL.TileSpacing = 'compact';
TL.Padding = 'compact';

% Create a vector that holds the info for which comparisons are
% statistically significantly different
rej_null_hyp = p_values_most_conservative_comps_scaled < 0.05;

% Write out this info to the command window so the user can quickly see
% which ones are statistically signifantly different
disp('The following are statistically significantly different:');
for i = find(rej_null_hyp)'
    disp([p_values_most_conservative_comps_sorted_names{i}]);
end