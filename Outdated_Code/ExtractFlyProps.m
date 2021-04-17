function ExtractFlyProps(inputFile) 

% Open video reader and figure out dimensions and frames of video
v = VideoReader(inputFile);
frameCount = floor(v.Duration * v.FrameRate) - 1;
resolution = [v.Width v.Height];

% Load or compute the rough trace of the individual cassettes
if isfile(strcat(erase(inputFile,'.avi'), '_cassettes.png'))
    cassetteImg = strcat(erase(inputFile,'.avi'), '_cassettes.png');
    grayCasImg = imread(cassetteImg);
    BWCasImg = imbinarize(grayCasImg);
    LabeledOrArray = bwlabel(BWCasImg);
else
    LabeledOrArray = CassetteTracer(inputFile);
end

% Compute the number of cassettes that were used throughout the experiment
numCassettes = max(LabeledOrArray, [], 'all');

% Give warning message if <40 cassettes are identified in an experiment
if numCassettes < 40
    imshow(LabeledOrArray);
    error = input('Warning: <40 cassettes used. Continue? y/n\n','s');
        if error == 'y'
            % Continue
        else
            input('Ctrl+C to cancel.\n');
        end
end

% Setup the structure to save the data
frame = struct('CassetteID', cell(1,frameCount));
sortedFrame = struct('CassetteID', cell(1,frameCount));

% Convert video to grayscale, binarize it, remove small blobs, find each
% fly, save each fly's size, centroid position, orientation, flyID, for
% each frame
for i = 1:frameCount
    gray = rgb2gray(read(v, i));
    BW = imbinarize(gray);
    BW2 = bwareaopen(BW,100);
    stats = regionprops('struct', BW2, 'Area', 'Centroid', 'Orientation');
    NumFliesVec = size(struct(stats));
    NumFlies(i) = NumFliesVec(1);
    % Make sure there aren't more flies detected than there are cassettes
    % If there are, set all values to 0 (would be nice to change this)
    if NumFlies(i) > numCassettes
        frame(i).CassetteID(1:numCassettes) = ...
            (struct('Area', {0},'Centroid', {[0 0]}, 'Orientation',{0})).';
    else
        frame(i).CassetteID = struct(stats);
    end
    % Make sure there are always labels for all flies, even during frames
    % where they are not being tracked
    if NumFlies(i) < numCassettes
        frame(i).CassetteID(end+1 : end+(numCassettes-NumFlies(i))) = ...
            struct('Area', {0},'Centroid', {[0 0]}, 'Orientation',{0});
    end
    % Add new fields to the structure (frameNum and CassetteNum)
    for j = 1:numCassettes
        frame(i).CassetteID(j).('frameNum') = i;
        % Determine CassetteNum by dotting centroid into labeled cassette
        % image
        if frame(i).CassetteID(j).Centroid(2) > 0
            tempCentroid = zeros(resolution(2), resolution(1));
            tempCentroid(floor(frame(i).CassetteID(j).Centroid(2)),...
                floor(frame(i).CassetteID(j).Centroid(1))) = 1;
            tempCentroid(floor(frame(i).CassetteID(j).Centroid(2)),...
                ceil(frame(i).CassetteID(j).Centroid(1))) = 1;
            tempCentroid(ceil(frame(i).CassetteID(j).Centroid(2)),...
                floor(frame(i).CassetteID(j).Centroid(1))) = 1;
            tempCentroid(ceil(frame(i).CassetteID(j).Centroid(2)),...
                ceil(frame(i).CassetteID(j).Centroid(1))) = 1;
            LabeledCentroid = LabeledOrArray.*tempCentroid;
            frame(i).CassetteID(j).('CassetteNum') = ...
                max(LabeledCentroid,[],'all');
        else
            frame(i).CassetteID(j).('CassetteNum') = 0;
        end
    end
    
    % Fill in CassetteID for the flies not tracked in a frame
    missingFlies = setdiff(1:numCassettes,[frame(i).CassetteID.CassetteNum]);
    for j = 1:(numCassettes-NumFlies(i))
        frame(i).CassetteID(NumFlies(i)+j).CassetteNum = missingFlies(j);
    end
    
    % Make a structure sorted by CassetteID
    assignin('base', 'frame', frame);
    ordCassettes = table2struct(sortrows(struct2table(frame(i).CassetteID),...
        'CassetteNum'));
    sortedFrame(i).CassetteID = struct(ordCassettes);
end

% Make final structure which now uses the CasseteID (equivalent to flyID
% at this step) as the main object and all other information as
% properties of each FlyID
FlyID = struct('FlyProps',cell(1,numCassettes));
for i = 1:numCassettes
    FlyID(i).FlyProps = struct('Centroid',{cell(frameCount,1)},...
        'Orientation',{cell(frameCount,1)},'Frame',{cell(frameCount,1)},...
        'Area',{cell(frameCount,1)},'CassetteID',{cell(frameCount,1)});
    for j = 1:frameCount
        FlyID(i).FlyProps.Centroid(j) = ...
            {sortedFrame(j).CassetteID(i).Centroid};
        FlyID(i).FlyProps.Orientation(j) = ...
            {sortedFrame(j).CassetteID(i).Orientation};
        FlyID(i).FlyProps.Frame(j) = ...
            {sortedFrame(j).CassetteID(i).frameNum};
        FlyID(i).FlyProps.Area(j) = ...
             {sortedFrame(j).CassetteID(i).Area};
        FlyID(i).FlyProps.CassetteID(j) = ...
             {sortedFrame(j).CassetteID(i).CassetteNum};
    end
end
    
% Save the structures to the workspace
assignin('base', 'sortedFrame', sortedFrame);
assignin('base', 'FlyID', FlyID);

end
