function MultiAnalyzeCrossStats(NumGenotypes)

% Close all figures to avoid any conflicts
close all

% Initialize vectors that hold the names of the files to be compared
WS_names = cell(NumGenotypes, 1);
WS_names_NoPath = cell(NumGenotypes, 1);

% Have the user select the genotype folders they want to compare
for genotypeCounter = 1:NumGenotypes
    folderPath = uigetdir('C:\Users\clarklab\Joe\Gap_Crossing\Data');
    % Fill in WS_names and WS_names_NoPath
    WS_names{genotypeCounter} = folderPath;
    WS_names_NoPath{genotypeCounter} = ...
        erase(folderPath,'C:\Users\clarklab\Joe\Gap_Crossing\Data\');
end

% Now go through and load in all the files
for genotypeCounter = 1:NumGenotypes
    % Load in the files
    WS = load([WS_names{genotypeCounter}, '\', ...
               WS_names_NoPath{genotypeCounter}, '_WS_Combined.mat']);
    % Remove the extra WS layer that comes from using load above
    WS = WS.WS;
    
    % Now run AnalyzeCrossStats
    WS = AnalyzeCrossStats(WS);
    
    % Clear WS for the next loop iteration
    clear WS
end

end