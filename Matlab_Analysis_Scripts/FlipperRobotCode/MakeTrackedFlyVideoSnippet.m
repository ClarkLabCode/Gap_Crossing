fileName = 'C:\Users\clark\Joe\Gap_Crossing\Data\All_Raw_Videos\2022-06-26 11-51-52_IsoD1_lighting_top_bottom_Exp1.mov';
vR = VideoReader(fileName);
vW = VideoWriter('Tracked_Flies_Video_for_Presentation.mov');
vW.FrameRate = 15;

numFlies = length(WS.FlipBinnedFlyStruct.ExpNum);
flipNum = 2;
snippetVec = 5:200;
for i_flies = 1:numFlies
    CentX(i_flies,:) = WS.FlipBinnedFlyStruct.ExpNum(i_flies).BehavData.OddFlips(flipNum).CentroidX(snippetVec);
    CentY(i_flies,:) = WS.FlipBinnedFlyStruct.ExpNum(i_flies).BehavData.OddFlips(flipNum).CentroidY(snippetVec);
    Frames = WS.FlipBinnedFlyStruct.ExpNum(i_flies).BehavData.OddFlips(flipNum).AbsoluteTime(snippetVec);
end

open(vW);

for i_frame = 1:length(snippetVec)
    img = read(vR,Frames(i_frame));
    img_withCirc = img;
    for i_flies = 1:numFlies
        img_withCirc = insertShape(img_withCirc,"circle",[CentX(i_flies,i_frame) CentY(i_flies,i_frame) 40],LineWidth=5,Color='red');
    end
    writeVideo(vW,img_withCirc);
end

close(vW);