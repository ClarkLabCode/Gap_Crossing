% Combines WS from several experiments into one common WS with a combined
% FlipBinnedFlyStruct that has the ExpNum layer removed

% Input:
% NumExps = Number of experiments to combine

function WS = CombineExpStructs(NumExp)

% Navigate to the data folder and let the user select the folder containing
% the data for the genotype of interest. This isn't necessary, but it saves
% the user time by helping automatically navigate folders later on.
flipperRobotCodePath = pwd;
dataPath = [flipperRobotCodePath,'\..\..\Data\'];
cd(dataPath);
dataPath = pwd;
dataPath = [dataPath, '\'];
cd(flipperRobotCodePath);
GenotypeDirectory = uigetdir(dataPath,'Select Genotype Folder');
cd(GenotypeDirectory);

% Initialize the combined WS and FBFS
WS = struct();
FBFS = [];

% Ask user to select all the experiments to be combined together
for expCounter = 1:NumExp
    
    % Request each WS file and set up the file path appropriately
    [file, path] = uigetfile('*.*', ['Select WS from Exp ', num2str(expCounter),...
        ' in the Fly_Structure folder'], ...
        [GenotypeDirectory, '\Experiment_', num2str(expCounter),'\5_Fully_Analyzed']);
    absFilePath = fullfile(path,file);
    inputFileName = absFilePath;
    
    % Load in WS into tempWS
    tempWS = load(inputFileName);
    
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

% Navigate back to the data directory
cd(flipperRobotCodePath);

end
