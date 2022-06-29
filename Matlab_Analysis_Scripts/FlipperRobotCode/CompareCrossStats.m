function CompareCrossStats(NumGenotypes)

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

% Now go through and load in all the figures
for genotypeCounter = 1:NumGenotypes
    % Load in the file
    openfig([WS_names{genotypeCounter},'\',WS_names_NoPath{genotypeCounter}, ...
          '_Proper_Crossing_Plot.fig']);
    openfig([WS_names{genotypeCounter},'\',WS_names_NoPath{genotypeCounter}, ...
          '_Crossing_Plot.fig']);
end

% Tile the figures on the screen (works well for 2-3 genotypes)
for genotypeCounter = 1:NumGenotypes
    movegui(figure(2*(genotypeCounter-1)+1), ...
            [10 + 600*(genotypeCounter-1) 560])
    movegui(figure(2*(genotypeCounter)), ...
            [10 + 600*(genotypeCounter-1) 50])
end
    
end