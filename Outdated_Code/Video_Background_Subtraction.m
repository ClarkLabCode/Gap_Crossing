%% Subtracts the background of a video and saves the new video
%% Must be used either with Video_Background_Finder.m or with its output
function Video_Background_Subtraction(inputFile, outputFile, nth_frame)

% Determine if background image exists or must be calculated
if isfile(strcat(erase(inputFile,'.mp4'), '_bg.png'))
    bg_img_name = strcat(erase(inputFile,'.mp4'), '_bg.png');
else
    bg_img_name = Video_Background_Finder(inputFile, nth_frame);
end

% Start read and write functions for new video
reader = VideoReader(inputFile);
writer = VideoWriter(outputFile);

% Match the properties of the output video to the input video
writer.FrameRate = reader.FrameRate;
open(writer);

% Copy each frame with the background subtracted out
while hasFrame(reader)
   img = imread(bg_img_name) - readFrame(reader);
   writeVideo(writer,img);
end

close(writer);

end
