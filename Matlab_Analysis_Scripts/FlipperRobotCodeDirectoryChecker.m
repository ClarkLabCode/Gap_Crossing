% Makes sure we are in the correct starting directory to use the
% FlipperRobot scripts by first trying to navigate the expected directory
% structure and then asking the user to locate the correct directory if the
% directory structure was not as expected. Then runs a verification that
% the user selected directory is indeed the expected directory. Once the
% directory is confirmed to be correct, it navigates to the directory and
% adds it to the path to enable the use of all scripts in the directory.

try
    % Check if we're in the correct directory by going up 3 directory
    % levels from the current directory and then navigating back in by the
    % expected 3 directory levels. Then add this to the path.
    cd ..\..\..\Gap_Crossing\Matlab_Analysis_Scripts\FlipperRobotCode
    addpath(pwd)
catch
    % If the above attempt to navigate up then down 3 directories didn't
    % work, ask the user to select the FlipperRobotCode directory
    flipperRobotCodeDirectoryPath = ...
        uigetdir('','Please select the FlipperRobotCode folder within \Gap_Crossing\Matlab_Analysis_Scripts\');
    % Check that the user's input directory is indeed correct by going up 3
    % directory levels from the selected directory and then navigating back
    % in by the expected 3 directory levels. If correct, add to path.
    try
        cd([flipperRobotCodeDirectoryPath, '\..\..\..\Gap_Crossing\Matlab_Analysis_Scripts\FlipperRobotCode']);
        addpath(flipperRobotCodeDirectoryPath);
    % If the selected directory was wrong, throw an error.
    catch
        error('The FlipperRobotCode folder was incorrectly chosen. Please try again.');
    end
end