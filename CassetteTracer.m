%% Generates a binary image that approximates the positions of the 
%% cassettes in a video by tracing the path that flies took during an
%% experiment. This binary image can be used to assign FlyIDs to each
%% fly in an experiment (use ExtractFlyProps.m for this).
function LabeledOrArray = CassetteTracer(inputFile)

% Open video reader and figure out dimensions and frames of video
v = VideoReader(inputFile);
frameCount = floor(v.Duration * v.FrameRate) - 1;
resolution = [v.Width v.Height];

% Make an image that will track all pixels that were on during video
orArray = zeros(resolution(2), resolution(1));

% Convert video to grayscale, binarize, remove small blobs, then fill in
% image of all pixels that were on during video (should roughly trace the
% cassettes)
for i = 1:frameCount
    gray = rgb2gray(read(v, i));
    BW = imbinarize(gray);
    BW2 = bwareaopen(BW,100);
    orArray = orArray|BW2;
end

% Give each traced cassette an ID number using bwlabel
LabeledOrArray = bwlabel(orArray)

% Show the image to make sure then save image to workspace and file
imshow(LabeledOrArray);
assignin('base','LabeledOrArray',LabeledOrArray);
imwrite(LabeledOrArray,strcat(erase(inputFile,'.avi'), '_cassettes.png'));

end