cd Training_Samples_Images\

rmdir('on_glass_crossing','s')
rmdir('off_glass_crossing','s')
% rmdir('lean_off_glass_crossing','s')
% rmdir('lean_on_glass_crossing','s')
% rmdir('ambiguous', 's')

mkdir('on_glass_crossing')
mkdir('off_glass_crossing')
% mkdir('lean_off_glass_crossing')
% mkdir('lean_on_glass_crossing')
mkdir('ambiguous', 's')

onCounter = 0;
offCounter = 0;
onLeanCounter = 0;
offLeanCounter = 0;
ambiguousCounter = 0;

for i = 1:length(LabeledCrossingVecLabel)
    if LabeledCrossingVecLabel(i) == 1
        cd off_glass_crossing\
        offCounter = offCounter + 1;
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(offCounter), '.png'])
        cd ..
%     else
    elseif LabeledCrossingVecLabel(i) == -1
        cd on_glass_crossing\
        onCounter = onCounter + 1;
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(onCounter), '.png'])
        cd ..
%     elseif LabeledCrossingVecLabel(i) == 0.5
%         cd lean_off_glass_crossing\
%         offLeanCounter = offLeanCounter + 1;
%         imwrite(LabeledCrossingMatInput(:,:,i),[num2str(offLeanCounter), '.png'])
%         cd ..
%     elseif LabeledCrossingVecLabel(i) == -0.5
%         cd lean_on_glass_crossing\
%         onLeanCounter = onLeanCounter + 1;
%         imwrite(LabeledCrossingMatInput(:,:,i),[num2str(onLeanCounter), '.png'])
%         cd ..
    else
        cd ambiguous
        ambiguousCounter = ambiguousCounter + 1;
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(ambiguousCounter), '.png'])
        cd ..
    end
end

cd ..