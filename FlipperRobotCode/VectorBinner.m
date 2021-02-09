% Takes a vector and bins the data by finding the average of every n elements
% Corresponds to Robot_Test20

function BinnedVec = VectorBinner(vec, n)

ElemCount = numel(vec);
% Removes the last data points in vec that don't fit into n elements
vecResh = reshape(vec(1:ElemCount - mod(ElemCount, n)), n, []);
BinnedVec = sum(vecResh)/n;

end

