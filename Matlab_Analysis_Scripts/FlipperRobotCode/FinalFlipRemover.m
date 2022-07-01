% Removes the row corresponding to the last flip for all flies because this
% row is always empty (because it is never a complete flip and incomplete
% flips are deemed empty).

% The final flip can sometimes give classifying errors after the data has  
% been processed through the neural net because the neural net loops
% through all frames of the video rather than just frames contained in
% completed flips. This was implemented this way because it drastically
% speeds up the rate at which the neural net can read in the video file, so
% it is worth adding in this extra error-prevention function for the sole
% purpose of keeping the code as efficient as possible.

function WS = FinalFlipRemover(WS)

% Port in the relevant fields from WS
FBFS = WS.FlipBinnedFlyStruct;
NumFlies = length(FBFS.ExpNum);

% Go through each fly and remove the row corresponding to the final flip
for flyCounter = 1:NumFlies
    FBFS.ExpNum(flyCounter).BehavData.OddFlips(end) = [];
    FBFS.ExpNum(flyCounter).BehavData.EvenFlips(end) = [];
end

% Update the fields in WS
WS.FlipBinnedFlyStruct = FBFS;

end
