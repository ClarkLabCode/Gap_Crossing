% LOADGCEXPWSNAMES Allows users to continue analysis from a checkpoint
%
%  This function is called in FullGapCrossingAnalysis to let users load in 
%  previous analysis that was saved in a checkpoint in FullGapCrossingAnalysis.
%  It asks the user to load in previously analyzed data one experiment at a
%  time and then outputs a message to the command window to confirm that
%  the file was loaded.

function WS_names = LoadGCExpWSNames(NumGCExp,LocalGCDirectoryPath)

% Initialize WS_names which holds the names for all the WS files
WS_names = cell(NumGCExp,1);

% Have the user select each file
for ExpCounter = 1:NumGCExp
    % Have the user select the file to load in and default to the data folder
    [file, path] = uigetfile('*.*','Select the WS struct file.',...
        [LocalGCDirectoryPath,'Data\']);
    pathMinusGCPath = erase(path,LocalGCDirectoryPath);
    WS_names{ExpCounter} = [pathMinusGCPath, file];
    % Output to command window a message that tells user which file was loaded
    disp(['Loaded ', path, file]);
end

end