% Combines WS from several experiments into one common WS with a combined
% FlipBinnedFlyStruct that has the ExpNum layer removed

% Inputs to varargin:
% 1) NumExps = Number of experiments to combine (can be numeric or 'All')
% 2) Genotype Folder = Directory to combine

function WS = CombineExpStructs(varargin)

% Determine how many inputs were given
if nargin == 1
    NumExp = varargin{1};
elseif nargin == 2
    NumExp = varargin{1};
    GenotypeDirectory = varargin{2};
else
    error('Unexpected number of inputs provided');
end

% Navigate to the data folder and let the user select the folder containing
% the data for the genotype of interest. This isn't necessary, but it saves
% the user time by helping automatically navigate folders later on.
flipperRobotCodePath = pwd;
dataPath = [flipperRobotCodePath,'\..\..\Data\'];
cd(dataPath);
dataPath = pwd;
dataPath = [dataPath, '\'];
cd(flipperRobotCodePath);
% If the user hasn't already passed the GenotypeDirectory to the function,
% pop up a GUI that has the user choose the Genotype Folder
if nargin ~= 2
    GenotypeDirectory = uigetdir(dataPath,...
        'Select Genotype Folder if combining one genotype. Otherwise select the Data folder.');
end
cd(GenotypeDirectory);

% Initialize the combined WS and FBFS
WS = struct();
FBFS = [];

% Check if NumExp was provided as 'All' or a number
if ~isnumeric(NumExp)
    % Set NumExp to number of experiments in folder
    if strcmpi('All',NumExp)
        subFolds = dir(GenotypeDirectory); % Grab everything in GenotypeDirectory
        subFolds = subFolds([subFolds.isdir]); % Only keep subfolders
        NumExp = size(subFolds,1)-2; % Subtract 2 because of '.' and '..' "directories"
        inputAll = 1; % Used later to bypass user GUIs when 'All' isn't the input
    % If something other than 'All' or a number was provided, throw error
    else
        cd(flipperRobotCodePath); % Move back to original directory before throwing error
        error('Expected first input to be an integer or ''All''');
    end
% Since NumExp was provided as a number, user did not request to auto
% input all experiments to be combined, so set inputAll to 0
else
    inputAll = 0;
end

% Initialize vector to hold all file names of structs being combined
inputFileNameVec = cell(NumExp,1);

% Ask user to select each experiment to be combined together
if inputAll == 0
    for expCounter = 1:NumExp
        % Request each WS file and set up the file path appropriately
        [file, path] = uigetfile('*.*', ['Select WS from Exp ', num2str(expCounter),...
            ' in the 5_Fully_Analyzed folder'], ...
            [GenotypeDirectory, '\Experiment_', num2str(expCounter),'\5_Fully_Analyzed']);
        % Combine everything and then feed it into inputFileNameVec
        absFilePath = fullfile(path,file);
        inputFileName = absFilePath;
        inputFileNameVec{expCounter} = inputFileName;
    end
% Auto-select all experiments to be combined together within GenotypeDirectory
else
    for expCounter = 1:NumExp
        % Grab everything inside Fully Analyzed in Experiment_expCounter in Genotype Folder
        fileStruct = dir([GenotypeDirectory, '\Experiment_', num2str(expCounter),'\5_Fully_Analyzed']);
        % Only keep the file and throw out the "directories" named '.' and '..'
        fileStruct = fileStruct(~[fileStruct.isdir]);
        % Pull out the name of the file
        file = fileStruct(1).name;
        % Pull out the path of the file
        path = [GenotypeDirectory, '\Experiment_', num2str(expCounter),'\5_Fully_Analyzed'];
        % Combine everything and then feed it into inputFileNameVec
        absFilePath = fullfile(path,file);
        inputFileName = absFilePath;
        inputFileNameVec{expCounter} = inputFileName;
    end
end
   
% Make a variable that toggles whether or not to prompt user after each
% warning as well as a variable that tracks whether or not there was a
% mismatch warning
skipContinuationPrompts = 0;
mismatchWarnings = 0;

% Go through and load in all the files
for expCounter = 1:NumExp    
    % Load in WS into tempWS
    tempWS = load(inputFileNameVec{expCounter});
    
    % Peel away the protective layer that comes from loading in a struct
    tempWS = tempWS.WS;
    
    % Fill in tempFBFS and remove the ExpNum layer
    tempFBFS = tempWS.FlipBinnedFlyStruct.ExpNum;
    % Also grab the skeleton coordinates of each fly's appropriate corridor
    for flyCounter = 1:size(tempFBFS,2)
        tempFBFS(flyCounter).IdData.Skeleton_x = tempWS.Skeleton_x(:,:,flyCounter);
        tempFBFS(flyCounter).IdData.Skeleton_y = tempWS.Skeleton_y(:,:,flyCounter);
    end
    
    % Concatenate tempFBFS into FBFS
    FBFS = [FBFS, tempFBFS];
    
    % Check that the experiments being combined aren't mismatching
    if expCounter ~= 1
        % Compare genotypes of experiments being combined
        if ~strcmpi(temp_genotype(find(~isspace(temp_genotype))), tempWS.genotype(find(~isspace(tempWS.genotype))))
            mismatchWarnings = 1;
            warning('Experiments with mismatching genotypes being combined');
            % If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
        % Compare directories of experiments being combined
        if ~strcmpi(temp_directoryName, tempWS.directoryName)
            mismatchWarnings = 1;
            warning('Experiments with mismatching directories being combined');% If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
        % Compare flip rates of experiments being combined
        if temp_flipRate            ~= tempWS.flipRate
            mismatchWarnings = 1;
            warning('Experiments with mismatching flip rates being combined');% If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
        % Compare experiment lengths of experiments being combined
        if temp_expLength            ~= tempWS.expLength
            mismatchWarnings = 1;
            warning('Experiments with mismatching lengths being combined');% If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
        % Compare temperatures of experiments being combined
        if temp_temperature         ~= tempWS.temperature
            mismatchWarnings = 1;
            warning('Experiments with mismatching termperatures being combined');% If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
        % Compare number of gaps of experiments being combined
        if temp_NumGaps             ~= tempWS.NumGaps
            mismatchWarnings = 1;
            warning('Experiments with mismatching number of gaps being combined');% If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
        % Compare neural net classification thresholds of experiments being combined
        if temp_threshProb          ~= tempWS.threshProb
            mismatchWarnings = 1;
            warning('Experiments with mismatching neural net classification settings being combined');% If the user has not requested to ignore continuation prompts,
            % ask them each time there is a mismatch if they want to
            % terminate the combining of the loaded experiments
            if skipContinuationPrompts ~= 1
                % Give user the option to terminate combining experiments
                if ~strcmpi('y',input('Continue? [y/n]\n', 's'))
                    error('User terminated combining experiments due to mismatch.')
                % If user chose to continue combining experiments, ask if they
                % want to no longer be prompted for continuing for subsequent
                % mismatches found in the combination process
                else
                    if strcmp('y',input('Ignore rest of mismatch warnings when combining these experiments? [y/n]\n','s'))
                        skipContinuationPrompts = 1;
                    end
                end
            end
        end
    end
    temp_genotype                    = tempWS.genotype;
    temp_directoryName               = tempWS.directoryName;
    temp_flipRate                    = tempWS.flipRate;
    temp_expLength                   = tempWS.expLength;
    temp_temperature                 = tempWS.temperature;
    temp_NumGaps                     = tempWS.NumGaps;
    temp_threshProb                  = tempWS.threshProb;
    
end

% Fill in the fields that are independent of whether or not there was
% mismatch between experiments being combined
WS.NumExp              = NumExp;
WS.NumFlies            = length(FBFS);
WS.FlipBinnedFlyStruct = FBFS;
WS.netOnOffAmbig       = tempWS.netOnOffAmbig;

% If there were no mismatches, fill in the fields that are common for all 
% experiments in WS
if ~mismatchWarnings
    WS.genotype            = tempWS.genotype;
    WS.directoryName       = tempWS.directoryName;
    WS.flipRate            = tempWS.flipRate;
    WS.expLength           = tempWS.expLength;
    WS.temperature         = tempWS.temperature;
    WS.NumGaps             = tempWS.NumGaps;
    WS.threshProb          = tempWS.threshProb;
% If there were mismatches, fill in the fields that may be different
% between experiments with "Potentially variable"
else
    WS.genotype            = 'Potentially variable';
    WS.directoryName       = 'Potentially variable';
    WS.flipRate            = 'Potentially variable';
    WS.expLength           = 'Potentially variable';
    WS.temperature         = 'Potentially variable';
    WS.NumGaps             = tempWS.NumGaps;
    WS.threshProb          = 'Potentially variable';
end
    
% Now save WS_Combined
% If there were no mismatches, then save it in the directoryName folder
if ~mismatchWarnings
    save([WS.directoryName,'_WS_Combined.mat'],'WS')
% If there were mismatches, ask the user where to save the combined WS and
% give them the choice of making a new folder for it
else
    if strcmpi('y',input('Do you want to save this combined WS into an already existing folder? [y/n]\n','s'))
        % Ask user to select the folder to save in if it's in an existing folder
        saveFolderLoc = uigetdir(dataPath,'Select folder to save in');
        cd(saveFolderLoc);
    else
        % Ask user the name of the new folder to make to hold the combined
        % WS in the data folder
        cd(dataPath)
        newFolderName = input('What should the name of the new folder be?\n','s');
        mkdir(newFolderName);
        cd(newFolderName);
    end

    % Ask user what to name the file in the chosen folder (new or pre-existing)
    combinedWSName = input('What should this combined WS be named in the selected folder?\n','s');
    save([combinedWSName,'_WS_Combined.mat'],'WS');
end

% Navigate back to the code directory
cd(flipperRobotCodePath);

end