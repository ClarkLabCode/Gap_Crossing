lengthInFlip = 240;
lengthBWFlip = 25;

expLength = 30*60*30;

numFlips = floor(expLength/(lengthInFlip+lengthBWFlip));

basisVec = zeros(lengthInFlip+lengthBWFlip,2);

basisVec(1:lengthInFlip,1) = -1;
basisVec(1:lengthInFlip,2) = 1;
basisVec((lengthInFlip+1):end,1) = (-1+2/(lengthBWFlip+1)):2/(lengthBWFlip+1):(1-2/(lengthBWFlip+1));
basisVec((lengthInFlip+1):end,2) = -1*((-1+2/(lengthBWFlip+1)):2/(lengthBWFlip+1):(1-2/(lengthBWFlip+1)));

basisVec1D = reshape(basisVec,[],1);

flipState = repmat(basisVec1D,[round(numFlips/2),1]);
flipState = flipState(1:expLength);

plot((1:expLength)/30,flipState,'k');
xlim([0,120]);
ylim([-1.2,1.2]);
xticks(0:30:120);
yticks([-3,3]);
ylabel('Flip Orientation');
xlabel('Time (s)');