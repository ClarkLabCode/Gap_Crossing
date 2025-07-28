% code to play with getting relative distances from angular velocities

% data = load('ParallaxSimsUpdated.mat');
% 
% %%
% 
% D1 = data.D1;
% D2 = data.D2;
% phiDot1 = data.phiDot1;
% phiDot2 = data.phiDot2;
% 
% %%
% 
% clear data;

%% find ratios and regions of the simulation

logDRatio = log10(D1./D2);
logPhiDotRatio = log10(abs(phiDot2./phiDot1)); 

% positive phiDot corresponds to back to front: simulating right eye with
% normal, right-handed coordinate system
s1 = find((phiDot1 < 0) .* (phiDot2 < 0)); % ftb | ftb
s2 = find((phiDot1 > 0) .* (phiDot2 > 0)); % btf | btf
s3 = find((phiDot1 < 0) .* (phiDot2 > 0)); % ftb | btf
s4 = find((phiDot1 > 0) .* (phiDot2 < 0)); % btf | ftb

%% sanity check on what's in the data...

disp(['ftb | ftb = ' num2str(length(s1)/length(D1))]);
disp(['btf | btf = ' num2str(length(s2)/length(D1))]);
disp(['ftb | btf = ' num2str(length(s3)/length(D1))]);
disp(['btf | ftb = ' num2str(length(s4)/length(D1))]);

%% some testing of the function performNonParametricFitting

xtest = randn(1,100000);
ytest = xtest + 1*randn(size(xtest));

performNonParametricFitting(xtest,ytest,500,0);
% note, this seems to work, fits well when noise is small, when noise is 1,
% then recovers 50% of the variance...

%% okay, now plot up with some real data

CD(1) = performNonParametricFitting(logPhiDotRatio(s1),logDRatio(s1),100,0);
CD(2) = performNonParametricFitting(logPhiDotRatio(s2),logDRatio(s2),100,0);
CD(3) = performNonParametricFitting(logPhiDotRatio(s3),logDRatio(s3),100,0);

% make a plot

figure; 
plot([1 2 3],CD,'k.','markersize',40);
ylabel('fraction variance in log(D2/D1) accounted for by log(phidot/phidot');
xlabel('motion type');
set(gca,'ylim',[0 1],'ytick',[0 0.5 1],'xtick',[1 2 3],'xlim',[0.5 3.5],...
    'xticklabel',{'FTB|FTB','BTF|BTF','FTB|BTF'},'xticklabelrotation',45);

%% mess with other options

% performNonParametricFitting(logPhiDotRatio,logDRatio,100,0);
% performNonParametricFitting(logPhiDotRatio(s3),10.^logDRatio(s3),100,0);
% performNonParametricFitting(10.^-logPhiDotRatio(s2),10.^logDRatio(s2),100,0);
% performNonParametricFitting(logPhiDotRatio(s3),10.^logDRatio(s3),100,0);





%% HELPER FUNCTIONS
function coeffDet = performNonParametricFitting(x,y,bins,pctlFlag, varargin)
    % this is a function to make a bunch of plots and analyses with x and y list as input, can eventually pass it a
    % list of handles for the different plots
    
    if pctlFlag
        Xedges = prctile(x,[0:100/bins:100]);
        Yedges = prctile(y,[0:100/bins:100]);
    else
        Xedges = linspace(min(x),max(x),bins+1);
        Yedges = linspace(min(y),max(y),bins+1);
    end
    Xcenters = Xedges(1:end-1)/2 + Xedges(2:end)/2; % average them
    Ycenters = Yedges(1:end-1)/2 + Yedges(2:end)/2; % average them

    if length(varargin)
        SILENT = varargin{1};
    else
        SILENT = 0;
    end

    Nxy = histcounts2(x,y,Xedges,Yedges)'; % transpose makes x values along columns in plots
    P_yGivenx = Nxy ./ repmat(sum(Nxy,1),[size(Nxy,1),1]);
    
    % start computing things usefully
    totVar = var(y);
    [Nx,dum,XbinIdx] = histcounts(x,Xedges);
    for ii=1:bins
        meanYgivenX(ii) = mean(y(XbinIdx == ii));
        yResidual(find(XbinIdx == ii)) = y(XbinIdx == ii) - meanYgivenX(ii);
        xMean(ii) = mean(x(XbinIdx == ii));
    end
    remainingVar = var(yResidual);
    coeffDet = (totVar - remainingVar)/totVar;
    disp(['yhat = E(y | x) has coefficient of determination of ' num2str(coeffDet)]);

    % redo the above calculation only for phi ratios in range 1/10 to 10
    ch = find(abs(x)<1); % run on this subset
    totVar = var(y(ch));
    [Nx,dum,XbinIdx] = histcounts(x(ch),Xedges);
    for ii=1:bins
        meanYgivenX(ii) = mean(y(ch(XbinIdx == ii)));
        yResidual(find(ch(XbinIdx == ii))) = y(ch(XbinIdx == ii)) - meanYgivenX(ii);
        xMean(ii) = mean(x(ch(XbinIdx == ii)));
    end
    remainingVar = var(yResidual(ch));
    coeffDet10 = (totVar - remainingVar)/totVar;
    disp(['within factor of 10: yhat = E(y | x) has coefficient of determination of ' num2str(coeffDet10)]);

    
    if ~SILENT
        % begin by plotting marginal distributions
        figure; 
        subplot(2,1,1);
        histogram(x,Xedges);
        xlabel('x values');
        ylabel('count');
        subplot(2,1,2);
        histogram(y,Yedges);
        xlabel('y values');
        ylabel('count');
        drawnow;
    
        % plot joint distribution
        figure;
        imagesc(log10(Nxy+1));
        xlabel('x value');
        ylabel('y value');
        colorbar;
        set(gca,'ydir','normal');
        % do x and y labels appropriately to get all integers in interval
        Xintvals = ceil(min(Xcenters)):floor(max(Xcenters));
        Yintvals = ceil(min(Ycenters)):floor(max(Ycenters));
        set(gca,'xtick',interp1(Xcenters,[1:length(Xcenters)],Xintvals),'XTickLabel',Xintvals);
        set(gca,'ytick',interp1(Ycenters,[1:length(Ycenters)],Yintvals),'YTickLabel',Yintvals);
        % colormap('copper')
        drawnow;
        
        % plot conditional distribution
        figure;
        imagesc(P_yGivenx);
        xlabel('x value');
        ylabel('y value');
        colorbar;
        set(gca,'ydir','normal');
        set(gca,'xtick',interp1(Xcenters,[1:length(Xcenters)],Xintvals),'XTickLabel',Xintvals);
        set(gca,'ytick',interp1(Ycenters,[1:length(Ycenters)],Yintvals),'YTickLabel',Xintvals);
        drawnow;
        
        % plot out E(y | x) vs. x
        figure;
        plot(xMean,meanYgivenX,'.-');
        xlabel('x value');
        ylabel('expected(y | x) in bin');
        title('plot of nonlinearity')
    end
    

end




