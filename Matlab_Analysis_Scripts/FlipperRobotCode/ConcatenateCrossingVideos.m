% Using this to generate video files that are a concatenation of all the
% crossing events used for training the neural net originally
% Useful for providing a training set for DLC limb labeling (Kate's project)

% Grab all the video clips
listOfVids = dir('D:\Old_Files\Gap_Crossing_Backup_01_19_22\Gap_Crossing\Data\All_Raw_Videos\Training_Videos');
listOfVids = listOfVids(3:end); % Filter out the . and .. directories

% Input number of different flies that were in training set
numFlies = 7;

% Initialize a vector that holds which video number in listOfVids is the
% last video corresponding to each fly
lastVidOfFly = zeros(numFlies,1);

% Go through the files and determine which video is the last that
% corresponds to each fly by sifting through the video file names
for flyCounter = 1:numFlies
    for vidCounter = 1:length(listOfVids)
        if contains(listOfVids(vidCounter).name,['Fly_',num2str(flyCounter)])
            lastVidOfFly(flyCounter) = vidCounter;
        end
    end
end

% Now just add 1 as the first entry of lastVidOfFly so it can easily be
% used for indexing
lastVidOfFly(end+1) = 1;
lastVidOfFly = circshift(lastVidOfFly,1);

% Now go through and open each video file corresponding to each fly and
% concatenate them all together into the new video file being written
% We write a separate concatenated video for each fly
for flyCounter = 1:numFlies
    % Initialize the video writer for each fly
    v_W = VideoWriter(['All_Crossings_Fly_',num2str(flyCounter)]);
    open(v_W);
    % Read all the videos corresponding to each fly and write to one video
    for vidCounter = lastVidOfFly(flyCounter):lastVidOfFly(flyCounter+1)
        v_R = VideoReader([listOfVids(vidCounter).folder, '\', listOfVids(vidCounter).name]);
        % Take all frames from each video and write them into the big one
        while hasFrame(v_R)
            writeVideo(v_W,readFrame(v_R));
        end
    end
    % Close the video writer after finishing each fly
    close(v_W);
end
