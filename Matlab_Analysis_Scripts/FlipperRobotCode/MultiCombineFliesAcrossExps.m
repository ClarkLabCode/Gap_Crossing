% Combines WS from all experiments of a given genotype into one common WS 
% with a combined FlipBinnedFlyStruct that has the ExpNum layer removed
% Can be run on several genotypes at a time

function MultiCombineFliesAcrossExps(varargin)

% Check to make sure that any inputs given to function make sense
if nargin == 1
    if isa(varargin{1}, double)
        NumGenotypes = varargin{1};
    else
        error('Expected input to be an integer.')
    end
elseif nargin > 1
    error('Expected a single input.')
end

% Move to the data folder as base directory and set the data path
cd 'C:\Users\clarklab\Joe\Gap_Crossing\Data\'
dataPath = 'C:\Users\clarklab\Joe\Gap_Crossing\Data\';

% Grab all the subdirectories in the folder by chceking the content in the
% folder and only grabbing content that is a directory and skipping the
% first two directories (which are always '.' and '..')
DataDir = dir(pwd);
DataDir = DataDir([DataDir.isdir]);
DataDir = DataDir(3:end);

% Now display all the subdirectory names and let the user select any
% subdirectories they want to combine experiments within
disp('Select the directories you want to combine experiments within.')
DesiredFolders = listdlg('ListString', {DataDir.name});

% Check to make sure that the number of folders selected matches the number
% of genotypes user input if user gave an input at function call
if nargin == 1
    if length(DesiredFolders) ~= NumGenotypes
        error('Wrong number of directories selected.')
    end
end

% Now convert the dialog input into an actual array of strings and append
% the absolute data path to all the inputs
genotypeFolderNames_all = {DataDir(DesiredFolders).name};
genotypeFolderNames_all = strcat(dataPath, genotypeFolderNames_all);

% Now cycle through each subdirectory and combine the experiments
for subdirCounter = 1:length(genotypeFolderNames_all)
    % Switch folders into the appropriate subdirectory for the genotype
    cd(genotypeFolderNames_all{subdirCounter});
    
    % Determine how many experiments were done for this genotype
    currDir = dir(pwd);
    NumExp = sum([currDir(~ismember({currDir.name},{'.','..'})).isdir]);

    % Initialize the combined WS and FBFS
    WS = struct();
    FBFS = [];
    
    % Go through all experiments within this genotype
    for expCounter = 1:NumExp
    
        % Grab the name of the fully analyzed file within each experiment
        fullyAnalyzedFilePath = [genotypeFolderNames_all{subdirCounter}, ...
            '\Experiment_', num2str(expCounter),'\5_Fully_Analyzed\'];
        fullyAnalyzedFileDir = dir(fullyAnalyzedFilePath);
        % The 3 below is to skip over the dummy directories '.' and '..'
        fullyAnalyzedFile = [fullyAnalyzedFilePath, fullyAnalyzedFileDir(3).name];

        % Load in WS into tempWS
        tempWS = load(fullyAnalyzedFile);

        % Peel away the protective layer that comes from loading in a struct
        tempWS = tempWS.WS;

        % Fill in tempFBFS and remove the ExpNum layer
        tempFBFS = tempWS.FlipBinnedFlyStruct.ExpNum;

        % Concatenate tempFBFS into FBFS
        FBFS = [FBFS, tempFBFS];

        % Check that the experiments being combined aren't mismatching
        if expCounter ~= 1
            if temp_genotype            ~= tempWS.genotype
                error('Mismatching experiments being combined');
            end
            if temp_directoryName       ~= tempWS.directoryName
                error('Mismatching experiments being combined');
            end
            if temp_flipRate            ~= tempWS.flipRate
                error('Mismatching experiments being combined');
            end
            if temp_temperature                ~= tempWS.temperature
                error('Mismatching experiments being combined');
            end
            if temp_NumGaps             ~= tempWS.NumGaps
                error('Mismatching experiments being combined');
            end
            if temp_threshProb          ~= tempWS.threshProb
                error('Mismatching experiments being combined');
            end
        end
        
        % Save the previous data to be able to compare
        temp_genotype                    = tempWS.genotype;
        temp_directoryName               = tempWS.directoryName;
        temp_flipRate                    = tempWS.flipRate;
        temp_temperature                 = tempWS.temperature;
        temp_NumGaps                     = tempWS.NumGaps;
        temp_threshProb                  = tempWS.threshProb;

    end

    % Fill in the fields that are common for all experiments in WS
    WS.genotype            = tempWS.genotype;
    WS.directoryName       = tempWS.directoryName;
    WS.flipRate            = tempWS.flipRate;
    WS.expLength           = tempWS.expLength;
    WS.temperature         = tempWS.temperature;
    WS.NumGaps             = tempWS.NumGaps;
    WS.NumExp              = NumExp;
    WS.NumFlies            = length(FBFS);
    WS.FlipBinnedFlyStruct = FBFS;
    WS.netOnOffAmbig       = tempWS.netOnOffAmbig;
    WS.threshProb          = tempWS.threshProb;

    % Now save WS_Combined
    save([WS.directoryName,'_WS_Combined.mat'],'WS')

end

end
