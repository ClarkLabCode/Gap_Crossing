% Function that lets user load in previous analysis into WS_all

function WS_all = LoadGCExpWS(NumGCExp)

% Initialize WS_all
WS_all = cell(NumGCExp,1);

% Initialize a variable that will hold all the file paths
absFilePath_all = cell(NumGCExp,1);

% Have the user select each file
for ExpCounter = 1:NumGCExp
    % Have the user select the file to load in and default to the data folder
    [file, path] = uigetfile('*.*','Select the WS struct file.','C:\Users\clarklab\Joe\Gap_Crossing\Data\');
    absFilePath = fullfile(path,file);
    absFilePath_all{ExpCounter} = absFilePath;
    % Output to command window a message that tells user which file was loaded
    disp(['Loaded ', absFilePath]);
end

% Load each file into the struct temp_WS then fill WS_all by unraveling temp_WS
for ExpCounter = 1:NumGCExp
    temp_WS = open(absFilePath_all{ExpCounter});
    WS_all{ExpCounter} = temp_WS.WS;
end

end