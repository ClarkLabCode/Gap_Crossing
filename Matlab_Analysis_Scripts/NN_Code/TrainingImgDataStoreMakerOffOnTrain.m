cd Off_On_Training_Samples_Images\

rmdir('on_glass_crossing','s')
rmdir('off_glass_crossing','s')

mkdir('on_glass_crossing')
mkdir('off_glass_crossing')

onCounter = 0;
offCounter = 0;
onLeanCounter = 0;
offLeanCounter = 0;
ambiguousCounter = 0;

for i = 1:round(length(LabeledCrossingVecLabel))
    if LabeledCrossingVecLabel(i) > 0
        cd off_glass_crossing\
%         offCounter = offCounter + 1;
%         imwrite(LabeledCrossingMatInput(:,:,i),[num2str(offCounter), '.png'])
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(i, '%05.0f'), '.png'])
        cd ..
    elseif LabeledCrossingVecLabel(i) < 0
        cd on_glass_crossing\
%         onCounter = onCounter + 1;
%         imwrite(LabeledCrossingMatInput(:,:,i),[num2str(onCounter), '.png'])
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(i, '%05.0f'), '.png'])
        cd ..
    end
end

cd ..