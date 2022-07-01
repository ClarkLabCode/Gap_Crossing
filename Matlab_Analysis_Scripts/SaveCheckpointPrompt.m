% Pops up a GUI that asks user to select whether to continue since they
% have reached a save checkpoint

function SaveCheckpointResponse = SaveCheckpointPrompt(contOpt, stopOpt)

% Create GUI to contain the prompt
stopFig = uifigure;

% Prompt user to choose between contOpt and stopOpt
SaveCheckpointResponse = uiconfirm(stopFig,...
    'Your analysis progress has been saved. Do you wish to continue?',...
    'Save Checkpoint', 'Options', {contOpt,stopOpt});

% Closer the GUI once response has been given
close(stopFig);

end