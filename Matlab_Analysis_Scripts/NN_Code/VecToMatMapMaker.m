% Reverse engineer which fly, what gap, and what event number each entry in
% LabeledCrossingVecLabel came from in CrossFrameIDLabels

% The correct answer gets output to VecToMatMap which is a cell array of
% the same length as LabeledCrossingVecLabel that holds as its entries the
% matrix entries in CrossFrameIDLabels that correspond to the element in
% LabeledCrossingVecLabel

VecCounter = 0;
VecToMatMap = cell(size(LabeledCrossingVecLabel));

for eventCounter = 1:19
    for gapCounter = 1:4
        for flyCounter = 1:7
            if CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) ~= 0
                VecCounter = VecCounter + 1;
                VecToMatMap{VecCounter} = [flyCounter, gapCounter, eventCounter];
            end
            if CrossFrameIDLabels(flyCounter,gapCounter+4,eventCounter) ~= 0
                VecCounter = VecCounter + 1;
                VecToMatMap{VecCounter} = [flyCounter, gapCounter+4, eventCounter];
            end
        end
    end
end

% Check to make sure that this actually worked. If it works, then
% errorCount will be equal to 0 after running this.
errorCount = 0;

for VecCounterCheck = 1:1000
    errorCount = errorCount + ...
        (LabeledCrossingVecLabel(VecCounterCheck) ~= ...
        CrossFrameIDLabels(...
            VecToMatMap{VecCounterCheck}(1),...
            VecToMatMap{VecCounterCheck}(2),...
            VecToMatMap{VecCounterCheck}(3)));
end