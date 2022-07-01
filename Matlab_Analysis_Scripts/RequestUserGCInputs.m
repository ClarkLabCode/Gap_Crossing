% Function that gathers necessary info for gap crossing (GC) experiments
% from user through assisted pop-up dialog boxes that try to parse
% available information from the file names of GC videos 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% | Assumes that the GC videos were recorded using OBS |
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function WS_all = RequestUserGCInputs(NumGCExp, StandardAnalysis)

% Initialize WS_all which will hold WS for all experiments
WS_all = cell(NumGCExp,1);

% Contains a list of all input parameters to be asked of user
ExpDlgList = {'Genotype','Save Directory Name','Input File Name','Experiment Repetition Number',...
    'Date Acquired','Time Acquired','Eclosion Date','Time Zone','Notes',...
    'Number of Corridors','Number of Gaps per Corridor','Gap Orientation',...
    'Flip Rate (sec)','Experiment Length (min)','Temperature (C)'};
% Contains a list of all non-standard analysis params to be asked of user
AnalysisDlgList = {'Neural Net Thresh (threshProb)','Fly Size Thresh (sizeThreshCutOff)', ...
    'Cassette Flip Buffer (indPosFrameBuffer)','Number of Pixels to Erode (erodePix)',...
    'Skeleton Example Coordinates', 'Skeleton Example Image',...
    'On/Off Glass Classification Neural Net'};

% Necessary title box for the input dialog boxes being made later
dlgtitle = 'Input';

% Set up dimensions of inputdlg box for experimental params
dimsExpH = [1,1,4,1,...
    1,1,1,1,5,...
    1,1,1,...
    1,1,1];
dimsExpW = ones(1,length(dimsExpH))*50;
dimsExp = [dimsExpH;dimsExpW]';

% Set up dimensions of inputdlg box for analysis params
dimsAnalysisH = [1,1,1,1,2,2,2];
dimsAnalysisW = ones(1,length(dimsAnalysisH))*50;
dimsAnalysis = [dimsAnalysisH;dimsAnalysisW]';

% Have the user select all video files
[fileAll, path] = uigetfile('*.*', 'Select Videos to Analyze', ...
                            'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\', ...
                            'MultiSelect', 'on');
absFilePathAll = fullfile(path,fileAll);
inputFileNameAll = absFilePathAll;

% Check to make sure that the appropriate number of videos were selected
if length(inputFileNameAll) ~= NumGCExp
    % When there's only one experiment being loaded in, the length of
    % inputFileNameAll is equal to the length of the file in characters, so
    % check to make sure that this isn't what's causing the issue
    if NumGCExp == 1
        % When only one experiment loaded in, inputFileName is a variable
        % of type char
        if ~ischar(inputFileNameAll)
            error('The wrong number of videos were selected. Please try again.');
        end
    % If there's a mismatch in video files selected and it's more than 1
    else
        error('The wrong number of videos were selected. Please try again.');
    end
end

% Give the dialog box pop ups for each video and populate fields
for vidCounter = 1:NumGCExp
    
    % Have to be careful if only loading in 1 experiment because of brace
    % indexing errors, so split into two cases
    if NumGCExp ~= 1
        % Load in the info we already have per experiment
        file = fileAll{vidCounter};
        inputFileName = inputFileNameAll{vidCounter};
    else
        % Load in the info we already have per experiment
        file = fileAll;
        inputFileName = inputFileNameAll;
    end
    
    % Check the file name for any info that we can grab
    if contains(file,'Exp')
        % If experiment number is marked in file name, grab it
        LocOfExpInName = strfind(file,'Exp');
        defExpNum = file(LocOfExpInName+3);
    else
        % If experiment number isn't marked in file name, default to 1
        LocOfExpInName = length(file)+2; % +2 is for consistency with defGenotype
        defExpNum = 1;
    end
    
    % Name format from OBS has first 10 digits being date of video
    defDateAcq = datetime(file(1:10));
    % Name format from OBS has time of video as 12th-19th digits
    defTimeAcq = file(12:19);
    
    % Assumes file name is '[date] [time]_[Directory Name]_Exp[ExpNum]_[Notes].mov'
    defDirectory = file(21:LocOfExpInName-2);
    defNotes = file(LocOfExpInName+5:end-4);
   
    % Define default inputs for inputdlg boxes
    defExpInput = {'',defDirectory,inputFileName,defExpNum,...
        datestr(defDateAcq),defTimeAcq,datestr(defDateAcq-1),'12',defNotes,...
        '7','4','Vertical',...
        '8','30','34'};
    defAnalysisInput = {'0.5','100','5','1',...
%         ['..\..\Values_Needed_For_All_Experiments\',...
%         'Compartment_Labeling_Example_Coords.mat'],...
%         ['..\..\Values_Needed_For_All_Experiments\',...
%         'Compartment_Mask_Labeler.png'],...
%         ['..\..\Values_Needed_For_All_Experiments\',...
%         'netOnOffAmbig.mat']};
        ['C:\Users\clarklab\Joe\Gap_Crossing\Values_Needed_For_All_Experiments\',...
        'Compartment_Labeling_Example_Coords.mat'],...
        ['C:\Users\clarklab\Joe\Gap_Crossing\Values_Needed_For_All_Experiments\',...
        'Compartment_Mask_Labeler.png'],...
        ['C:\Users\clarklab\Joe\Gap_Crossing\Values_Needed_For_All_Experiments\',...
        'netOnOffAmbig.mat']};

    % Pop up a dialog box to fill out all experimental info into the fields
    ExpParams = inputdlg(ExpDlgList,dlgtitle,dimsExp,defExpInput,'on');
    
    % If non-standard analysis is being run, pop up a dialog box for it
    if ~StandardAnalysis
        AnalysisParams = inputdlg(AnalysisDlgList,dlgtitle,dimsAnalysis,defAnalysisInput,'on');
    % If standard analysis is being run, just auto-fill the values
    else
        AnalysisParams = defAnalysisInput;
    end
    
    % Load in the neural net for classifying crossings
    tempNNStruct = load(AnalysisParams{7});
    FieldsInNNStruct = fieldnames(tempNNStruct);
    % Check that only one object was loaded
    if length(FieldsInNNStruct) ~= 1
        error('Unexpected result when loading neural net.');
    else
        netOnOffAmbig = tempNNStruct.(FieldsInNNStruct{1});
    end

    % Checks to make sure that no critical fields are empty
    if (isempty(ExpParams{1}) || isempty(ExpParams{2}) || isempty(ExpParams{3}))
        error('Genotype, Directory Name, and/or File Name field(s) empty.');
    end
    
    % Initialize an empty WS struct to fill in all the params for current experiment
    WS = struct();
    
    % Convert inputs into variables and save them in WS
    WS.genotype                         = ExpParams{1};
    WS.directoryName                    = ExpParams{2};
    WS.inputFileName                    = ExpParams{3};
    WS.ExpRepNum                        = str2double(ExpParams{4});
    WS.dateAcq                          = ExpParams{5};
    WS.timeAcq                          = ExpParams{6};
    WS.eclosionDate                     = ExpParams{7};
    WS.timeZone                         = ExpParams{8};
    WS.notes                            = ExpParams{9};
    WS.NumCorridors                     = str2double(ExpParams{10});
    WS.NumGaps                          = str2double(ExpParams{11});
    WS.GapOrientation                   = ExpParams{12};
    WS.flipRate                         = str2double(ExpParams{13});
    WS.expLength                        = str2double(ExpParams{14});
    WS.temperature                      = str2double(ExpParams{15});
    WS.threshProb                       = str2double(AnalysisParams{1});
    WS.sizeThreshCutOff                 = str2double(AnalysisParams{2});
    WS.indPosFrameBuffer                = str2double(AnalysisParams{3});
    WS.erodePix                         = str2double(AnalysisParams{4});
    WS.ExampleSkeletonCoordsFilePath    = AnalysisParams{5};
    WS.ExampleSkeletonImgFilePath       = AnalysisParams{6};
    WS.NeuralNetFilePath                = AnalysisParams{7};
    WS.netOnOffAmbig                    = netOnOffAmbig;
    
    % Now fill in the appropriate entry of WS_all with WS then clear WS
    WS_all{vidCounter} = WS;
    clear WS
    
end

end