cd Ambig_Unambig_Training_Samples_Images\

rmdir('ambig','s')
rmdir('unambig','s')

mkdir('ambig')
mkdir('unambig')

ambiguousCounter = 0;
unambiguousCounter = 0;

for i = 1:round(length(LabeledCrossingVecLabel))
    if abs(LabeledCrossingVecLabel(i)) == 1
        cd unambig\
%         unambiguousCounter = unambiguousCounter + 1;
%         imwrite(LabeledCrossingMatInput(:,:,i),[num2str(unambiguousCounter), '.png'])
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(i, '%05.0f'), '.png'])
        cd ..
    elseif abs(LabeledCrossingVecLabel(i)) == 0.5
        cd ambig\
%         ambiguousCounter = ambiguousCounter + 1;
%         imwrite(LabeledCrossingMatInput(:,:,i),[num2str(ambiguousCounter), '.png'])
        imwrite(LabeledCrossingMatInput(:,:,i),[num2str(i, '%05.0f'), '.png'])
        cd ..
    end
end

cd ..