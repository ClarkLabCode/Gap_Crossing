for setCounter = 1:4
for flyCounterInSet = 1:7
% setCounter
% flyCounterInSet
numFramesInFlipOdd = length(LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.OddFlips(LabelOrder(flyCounterInSet,setCounter)).Area);
numFramesInFlipEven = length(LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.EvenFlips(LabelOrder(flyCounterInSet+NumFlies,setCounter)).Area);
% length(LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.OddFlips(LabelOrder(flyCounterInSet,setCounter)).OffGlassYN);
LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.OddFlips(LabelOrder(flyCounterInSet,setCounter)).OffGlassYN = ...
[LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.OddFlips(LabelOrder(flyCounterInSet,setCounter)).OffGlassYN((end-numFramesInFlipOdd+1):end)];
% length(LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.EvenFlips(LabelOrder(flyCounterInSet+NumFlies,setCounter)).OffGlassYN);
LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.EvenFlips(LabelOrder(flyCounterInSet+NumFlies,setCounter)).OffGlassYN = ...
[LabeledFlyStruct.ExpNum(flyCounterInSet).BehavData.EvenFlips(LabelOrder(flyCounterInSet+NumFlies,setCounter)).OffGlassYN((end-numFramesInFlipEven+1):end)];
end
end