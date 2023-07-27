function [net, YPred, YValidation, imdsValidation, imdsTest] = TrainGlassNNClassifier(folderName, classNames, biasWeight, numFramesGrabbed)

labelDatasetPath = ['C:\Users\clarklab\Joe\Gap_Crossing\Data\All_Raw_Videos\', folderName];

imds = imageDatastore(labelDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

labelCount = countEachLabel(imds);
[numLabels, ~] = size(labelCount);

img = readimage(imds,1);
[y,x] = size(img);

fracTrainFiles = 0.8;
[imdsTrain,imdsValidation,imdsTest] = ...
    splitEachLabel(imds,fracTrainFiles^2, fracTrainFiles-fracTrainFiles^2, (1-fracTrainFiles));
%         ,'exclude', 'all');

weightOfBias = biasWeight;

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
    classificationLayer('classes', classNames, 'ClassWeights', [1-weightOfBias weightOfBias])];
        % ('classes', ["ambig" "unambig"], 'ClassWeights', [1-weightOfBias weightOfBias])


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