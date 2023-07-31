% COMBINEANDANALYZEALLEXPSCRIPTS Automatically analyzes all folders in selected directory
%
%  Upon running this script, the user will be prompted to select the data
%  folder to be fully analyzed. If there are any folders you do not want
%  analyzed, then edit the line(s) above the tic command to exclude folders
%  you want to skip.
%
%  This full analysis consists of combining all experiments within each
%  folder (CombineExpStructs), aligning all fly data within these combined 
%  experiments (AlignCorridors), then running crossing statistics analysis 
%  and saving the plots and data (AnalyzeCrossStats).

% Ask user to select the data folder
AllAnalyzedExpFolders = uigetdir('D:\', 'Choose Data folder within Gap_Crossing folder');
% Extract all folders from inside the Data folder
AllAnalyzedExpFolders = AllAnalyzedExpFolders([AllAnalyzedExpFolders.isdir]);
% Remove the folders that you don't want to analyze
AllAnalyzedExpFolders(1:4) = []; % Remove ., .., All Raw Videos, All Plots
AllAnalyzedExpFolders(15) = []; % Remove IsoD1 8 all 2 mm
AllAnalyzedExpFolders(16) = []; % Remove IsoD1 all 2 mm

% Use tic and toc to track how long it took to complete all analysis
tic
% Now go through and combine all the experiments selected and analyze them
for i = 1:size(AllAnalyzedExpFolders,1)
    % Grab all the necessary info to feed into the combine and analysis functions
    GenotypeFolder = AllAnalyzedExpFolders(i).name;
    GenotypeFolderPath = [AllAnalyzedExpFolders(i).folder,'\'];
    GenotypeFolderFullPath = fullfile(GenotypeFolderPath,GenotypeFolder);
    % Combine all experiments of the current genotype
    WS = CombineExpStructs('All',GenotypeFolderFullPath);
    % Then align all the data from the different experiments for a given genotype
    WS = AlignCorridors(WS);
    % Produce the cross stats and figure for each genotype
    WS = AnalyzeCrossStats(WS);
    % Close the figure you generated
    close all
    % Give a progress update to the user
    disp(['Completed ', GenotypeFolder, ', genotype #', num2str(i)]);
end
toc