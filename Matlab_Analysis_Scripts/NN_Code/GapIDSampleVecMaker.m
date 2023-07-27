imdsSampleVec = imdsTest_OnOff.Files;
for i = 1:200
    imdsSampleVec{i} = erase(imdsSampleVec{i},...
        'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Off_On_Training_Samples_Images\off_glass_crossing\0');
    imdsSampleVec{i} = erase(imdsSampleVec{i},...
        'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Off_On_Training_Samples_Images\on_glass_crossing\0');
    imdsSampleVec{i} = erase(imdsSampleVec{i},...
        '.png');
    imdsSampleVec{i} = str2num(imdsSampleVec{i});
end
imdsSampleVec = cell2mat(imdsSampleVec);

gapIDSampleVec = zeros(200,1);
for i = 1:200
    dummyVec = cell2mat(VecToMatMap(imdsSampleVec(i)));
    gapIDSampleVec(i) = dummyVec(2);
end
gapIDSampleVec = gapIDSampleVec - 4*(gapIDSampleVec > 4);