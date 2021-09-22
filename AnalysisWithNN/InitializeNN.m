labelDatasetPath = 'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Training_Samples_Images';
% labelDatasetPath = 'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Off_On_Training_Samples_Images';
% labelDatasetPath = 'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Ambig_Unambig_Training_Samples_Images';
% labelDatasetPath = 'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Off_On_Testing_Samples_Images';
% labelDatasetPath = 'C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\Ambig_Unambig_Testing_Samples_Images';
imds = imageDatastore(labelDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

% figure;
% perm = randperm(500,20);
% for i = 1:20
%     subplot(4,5,i);
%     imshow(imds.Files{perm(i)});
% end

labelCount = countEachLabel(imds);
[numLabels, ~] = size(labelCount);

img = readimage(imds,1);
[y,x] = size(img);

fracTrainFiles = 0.8;
[imdsTrain,imdsValidation] = splitEachLabel(imds,fracTrainFiles);
%             ,'randomize'

for weightOfBias = 1

layers = [
    imageInputLayer([y x 1])
    
    convolution2dLayer([y/numFramesGrabbed x], 10, 'stride', [y/numFramesGrabbed 1],'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(numFramesGrabbed, 'stride', numFramesGrabbed, 'Padding','same')
    
    convolution2dLayer([y/numFramesGrabbed round(x/25)], 10, 'stride', [y/numFramesGrabbed 1],'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(numFramesGrabbed, 'stride', numFramesGrabbed, 'Padding','same')
        
    convolution2dLayer([y 1],10, 'stride', [y 1], 'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(numLabels)
    softmaxLayer
    classificationLayer];
%         ('classes', ["ambig" "unambig"], 'ClassWeights', [1-weightOfBias weightOfBias])

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',75, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(imdsTrain,layers,options);

YPred = classify(net,imdsValidation);
YValidation = imdsValidation.Labels;

accuracy = sum(YPred == YValidation)/numel(YValidation);

weightOfBias
numDataUnambig = sum((YPred == 'unambig'));
fracDataUnambig = numDataUnambig/length(YPred)
numFalseUnambig = sum((YValidation == 'ambig') & (YPred == 'unambig'));
fracFalseUnambig = numFalseUnambig/numDataUnambig

end