% FINALFLIPREMOVER Removes the row corresponding to the last flip of an experiment
%
%  Removes the row corresponding to the last flip for all flies in order to
%  make it safely compatible with the NN's classification scheme.
% 
%  The final flip can sometimes give classifying errors after the data has  
%  been processed through the neural net because the neural net loops
%  through all frames of the video rather than just frames contained in
%  completed flips. This was implemented this way because it drastically
%  speeds up the rate at which the neural net can read in the video file, so
%  it is worth adding in this extra error-prevention function for the sole
%  purpose of keeping the code as efficient as possible. This error
%  prevention results in a net loss of ~0.5% of our overall data.

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
