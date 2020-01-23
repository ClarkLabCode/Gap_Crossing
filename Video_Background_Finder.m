%% Computes the background in a video and saves as a png
%% Does not work when video has more than 16 million frames
%% Averages every nth_frame
function bg_img_name = Video_Background_Finder(inputFile, nth_frame)

% Read in the video
v = VideoReader(inputFile);

% Set up arrays for the background image 
frame = floor(v.Duration * v.FrameRate) - 1;
resolution = [v.Width v.Height];
bg_array32 = zeros(resolution(2), resolution(1), 'uint32'); %32-bit array that holds the sum of all frames
bg_step = 0;

% Go through every 300th frame of the video and add it to bg_array32
while (bg_step+1)*nth_frame < frame
    bg_step = bg_step + 1;
    bg_frame = rgb2gray(read(v, bg_step*nth_frame));
    bg_frame32 = uint32(bg_frame);
    bg_array32 = bg_frame32 + bg_array32;
end

% Find average frame by dividing by number of summed frames
bg_array8 = uint8(bg_array32/bg_step);

% Convert average frame to RGB image (necessary in RGB for other functions)
bg_array(:,:,1) = bg_array8;
bg_array(:,:,2) = bg_array8;
bg_array(:,:,3) = bg_array8;

% Save background image as png
bg_img_name = strcat(erase(inputFile,'.mp4'), '_bg.png')
imwrite(bg_array, bg_img_name);
    
end