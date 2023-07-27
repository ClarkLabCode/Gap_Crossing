% Get dimensions of input matrix
[y, x, numSamples] = size(LabeledCrossingMatInput);

% Reshape each sample into vector rather than image
LabeledCrossingVecInput = reshape(LabeledCrossingMatInput, [], numSamples);

% Perform PCA
[coeff,score,latent,tsquared,explained,mu] = pca(LabeledCrossingVecInput, 'NumComponents', 1000);

% Regenerate data from PCA
LabeledCrossingVecInputPCA = score*coeff';

% Reshape back to image from vector
LabeledCrossingMatInputPCA = reshape(LabeledCrossingVecInputPCA,y,x,numSamples);