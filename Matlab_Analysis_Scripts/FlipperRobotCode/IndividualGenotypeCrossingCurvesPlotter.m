% INDIVIDUALGENOTYPECROSSINGCURVESPLOTTER Script to plot individual crossing curves
%
%  This script allows the user to plot individual crossing curves (typically 
%  of one particular genotype, but not necessarily). In order to successfully 
%  run this script, the cell arrays AllCrossStats, AllCrossStatsNames, and
%  AllCrossStatsGenotypes must all be loaded into the workspace.
%
%  Users are given the choice between standard plots (i.e., comparing like
%  genotypes to controls) or making any other comparison they choose. They
%  are also given the option to plot the synthetic control, both parental
%  controls, or all three controls when making standard comparison plots.

% Check that the user has loaded in AllCrossStats
if ~exist('AllCrossStats','var')
    error('Please load in AllCrossStats, AllCrossStatsNames, and AllCrossStatsGenotypes before running this script.')
end

% Ask user what kind of comparison they want to plot
genotypesCompared = ...
    questdlg('What genotypes do you want to compare?','Genotype Comparison',...
             '[Neuron] > +;shts;shR to control(s)','Other comparison(s)',...
             '[Neuron] > +;shts;shR to control(s)');

% Present different prompts based on which comparison type was selected
switch genotypesCompared
    % If doing the standard comparison with shts (2&3)
    case '[Neuron] > +;shts;shR to control(s)'
        % Ask user which genotype they want to plot from AllCrossStats
        genotype = input('Which genotype do you want to plot? Open AllCrossStatsNames and respond with the row number of the desired genotype.\n');
        % Check that the user input a proper response
        if ~isnumeric(genotype)
            error('Expected a number to the previous prompt');
        end

        % Use split empty > shts as control for split Gal4s
        if contains(AllCrossStatsNames{genotype,2}, 'split')
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
        % 1) Plot all controls
        % 2) Plot the synthetic control curve only
        % 3) Plot the two genetic control curves
        whichControlsToPlot = questdlg('What control curve(s) would you like plotted?', ...
            'Choice of plotted control(s)', ...
            'All controls','Synthetic','Both genetic controls','All controls');

        % Open a new figure
        figure

        % Plot the selected choice of controls
        switch whichControlsToPlot
            % Plot the synthetic control and both genetic controls and legend
            case 'All controls'
                hold on
                % Silencing experiment (red)
                shadedErrorBar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd,...
                     'lineprops', {'Color','r','Marker','none'});
                % sp empty (DBD) > shts control (dark gray)
                shadedErrorBar(1:0.5:2.5,...
                     mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd,...
                     'lineprops', {'Color',[0.7,0.7,0.7],'Marker','none'});
                % Gal4 > + control (light gray)
                shadedErrorBar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd,...
                     'lineprops', {'Color',[0.3,0.3,0.3],'Marker','none'});
                % Synthetic control (black)
                shadedErrorBar(1:0.5:2.5,...
                     synth_mean,...
                     synth_sem,...
                     'lineprops', {'Color','k','Marker','none'});
                legend({'Silenced','sp empty (DBD) > shts','Gal4 > +','Synthetic control'})
    
            % If synthetic chosen, then plot the synth curve and legend
            case 'Synthetic'
                hold on
                % Silencing experiment (red)
                shadedErrorBar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd,...
                     'lineprops', {'Color','r','Marker','none'});
                % Synthetic control
                shadedErrorBar(1:0.5:2.5,...
                    synth_mean,...
                    synth_sem,...
                    'lineprops', {'Color','k','Marker','none'});
                legend({'Silenced','Synthetic Control'});
    
            % If both genetic controls chosen, then plot both and legend
            case 'Both genetic controls'
                hold on
                % Silencing experiment (red)
                shadedErrorBar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,2}.AllUpVecProperCrossOverAllGapEventStd,...
                     'lineprops', {'Color','r','Marker','none'});
                % Gal4 > + control (gray)
                shadedErrorBar(1:0.5:2.5,...
                     mean(AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventRate),...
                     AllCrossStats{genotype,3}.AllUpVecProperCrossOverAllGapEventStd,...
                     'lineprops', {'Color',[.7 .7 .7],'Marker','none'});
                % sp empty (DBD) > shts control (black)
                shadedErrorBar(1:0.5:2.5,...
                    mean(AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventRate),...
                    AllCrossStats{controlGenotype,1}.AllUpVecProperCrossOverAllGapEventStd,...
                    'lineprops', {'Color','k','Marker','none'});
                legend({'Silenced','Gal4 > +','sp empty (DBD) > shts',})
        end
    
        hold off
    
        % Just adjust the axes and add the appropriate plot and axes titles now
        ylim([0,1]);
        xlim([0.8,2.7]);
        xticks(1:0.5:2.5);
        title(AllCrossStatsNames{genotype,2}(1:(strfind(AllCrossStatsNames{genotype,2}, '>')-2)));
        xlabel('Gap Width (mm)');
        ylabel('Proper Crossings / All Gap Events');

    % If comparing experiments other than the standard shts (2&3) comparisons
    case 'Other comparison(s)'
        % Ask the user to select the desired curves to compare
        Gal4UASIdx = input('Select the experimental line (typically Gal4 > UAS) by typing the ordered pair of that genotype in AllCrossStatsNames.\n');
        EmptyUASIdx = input('Select the first parental control line (typically sp empty (DBD) > UAS) by typing the ordered pair of that genotype in AllCrossStatsNames.\n');
        Gal4PlusIdx = input('Select the second parental control line (typically Gal4 > +) by typing the ordered pair of that genotype in AllCrossStatsNames.\n');
        IsoD1Idx = input('Select the wild type line (typically IsoD1) by typing the ordered pair of that genotype in AllCrossStatsNames.\n');

        % Construct the synthetic controls (this is a model where the effect of
        % each construct [UAS and Gal4] have an individual additive effect)
        % The mean of the synthetic control algebraically simplifies to <Gal4/+> + <UAS/+> - <+> 
        synth_mean = ...
            mean(AllCrossStats{Gal4PlusIdx(1),Gal4PlusIdx(2)}.AllUpVecProperCrossOverAllGapEventRate) + ...
            mean(AllCrossStats{EmptyUASIdx(1),EmptyUASIdx(2)}.AllUpVecProperCrossOverAllGapEventRate) - ...
            mean(AllCrossStats{IsoD1Idx(1),IsoD1Idx(2)}.AllUpVecProperCrossOverAllGapEventRate);
        % Ensure that the crossing probability is non-negative
        synth_mean = synth_mean.*(synth_mean>=0);
        % The error of the synthetic control is computed via error propogation
        synth_sem = sqrt(...
            (AllCrossStats{Gal4PlusIdx(1),Gal4PlusIdx(2)}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
            (AllCrossStats{EmptyUASIdx(1),EmptyUASIdx(2)}.AllUpVecProperCrossOverAllGapEventStd).^2 + ...
            (AllCrossStats{IsoD1Idx(1),IsoD1Idx(2)}.AllUpVecProperCrossOverAllGapEventStd).^2);

        % Open a new figure
        figure

        hold on
        % Experimental line (usually Gal4 > UAS) (red)
        shadedErrorBar(1:0.5:2.5,...
             mean(AllCrossStats{Gal4UASIdx(1),Gal4UASIdx(2)}.AllUpVecProperCrossOverAllGapEventRate),...
             AllCrossStats{Gal4UASIdx(1),Gal4UASIdx(2)}.AllUpVecProperCrossOverAllGapEventStd,...
             'lineprops', {'Color','r','Marker','none'});
        % Parental line 1 (usually sp empty (DBD) > shts) (dark gray)
        shadedErrorBar(1:0.5:2.5,...
             mean(AllCrossStats{EmptyUASIdx(1),EmptyUASIdx(2)}.AllUpVecProperCrossOverAllGapEventRate),...
             AllCrossStats{EmptyUASIdx(1),EmptyUASIdx(2)}.AllUpVecProperCrossOverAllGapEventStd,...
             'lineprops', {'Color',[0.7,0.7,0.7],'Marker','none'});
        % Parental line 2 (usually Gal4 > + control) (light gray)
        shadedErrorBar(1:0.5:2.5,...
             mean(AllCrossStats{Gal4PlusIdx(1),Gal4PlusIdx(2)}.AllUpVecProperCrossOverAllGapEventRate),...
             AllCrossStats{Gal4PlusIdx(1),Gal4PlusIdx(2)}.AllUpVecProperCrossOverAllGapEventStd,...
             'lineprops', {'Color',[0.3,0.3,0.3],'Marker','none'});
        % Synthetic control (black)
        shadedErrorBar(1:0.5:2.5,...
             synth_mean,...
             synth_sem,...
             'lineprops', {'Color','k','Marker','none'});
        % Populate the legend using the names of the experiments used
        legend({AllCrossStatsNames{Gal4UASIdx(1),Gal4UASIdx(2)},...
                AllCrossStatsNames{EmptyUASIdx(1),EmptyUASIdx(2)},...
                AllCrossStatsNames{Gal4PlusIdx(1),Gal4PlusIdx(2)},...
                'Synthetic control'})
        hold off
    
        % Just adjust the axes and add the appropriate plot and axes titles now
        ylim([0,1]);
        xlim([0.8,2.7]);
        xticks(1:0.5:2.5);
        title(AllCrossStatsNames{Gal4UASIdx(1),Gal4UASIdx(2)}(1:(strfind(AllCrossStatsNames{Gal4UASIdx(1),Gal4UASIdx(2)}, '>')-2)));
        xlabel('Gap Width (mm)');
        ylabel('Proper Crossings / All Gap Events');

end