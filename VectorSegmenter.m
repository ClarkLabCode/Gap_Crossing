% Segments a vector into numSegs segments and puts them into a matrix
% Row 1 of matrix is the first n elements, row 2 is second n elements, etc
% Corresponds to Robot_Test21

function SegVec = VectorSegmenter(vec, numSegs)

ElemCount = numel(vec);
SegLength = floor(ElemCount/numSegs);
% Removes the last data points in vec that don't fit into n elements
SegVec = reshape(vec(1:ElemCount - mod(ElemCount, SegLength)), SegLength, []);
SegVec = SegVec';

end