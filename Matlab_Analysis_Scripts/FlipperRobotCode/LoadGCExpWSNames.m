% Function that lets user load in previous analysis into WS_all

function WS_names = LoadGCExpWSNames(NumGCExp,LocalGCDirectoryPath)

% Initialize WS_names which holds the names for all the WS files
WS_names = cell(NumGCExp,1);

% Have the user select each file
for ExpCounter = 1:NumGCExp
    % Have the user select the file to load in and default to the data folder
    [file, ~] = uigetfile('*.*','Select the WS struct file.',...
        [LocalGCDirectoryPath,'Data\']);
    WS_names{ExpCounter} = file;
    % Output to command window a message that tells user which file was loaded
    disp(['Loaded ', LocalGCDirectoryPath, 'Data\All_Raw_Videos\', file]);
end

end