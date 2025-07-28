% Plot Y vs X in a histogram for all flies
for oddEven = 1:2;
figure
numBins = 400;
histX = reshape([IsoD1MatrixX(:,:,:,oddEven)],[],1);
histX = histX(~isnan(histX));
histY = reshape([IsoD1MatrixY(:,:,:,oddEven)],[],1);
histY = histY(~isnan(histY));
XBinLim = 10;
YBinLim = 25;
histogram2(histX,histY,-XBinLim:(XBinLim/numBins):XBinLim,-YBinLim:(2*YBinLim/numBins):YBinLim,'DisplayStyle','tile','LineStyle','none','ShowEmptyBins','on');
set(gca,'ColorScale','log');
% colorbar;
grid off
% Example of how to overlay the skeleton in the plots
xl = xlim;
yl = ylim;
hold on
horizShift = 0.15;
% a = plot([0;StraightenedSkelX(end/2+1:end,1);0],...
%          [StraightenedSkelY(19,1);StraightenedSkelY(end/2+1:end,1);StraightenedSkelY(36,1)],'k');
a = plot([StraightenedSkelX(:,oddEven);StraightenedSkelX(1,oddEven)]-horizShift,...
         [StraightenedSkelY(:,oddEven);StraightenedSkelY(1,oddEven)],'k');
uistack(a,'top');
uistack(a,'top');
hold off
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
cmap = colormap('gray');
colormap(flip(cmap,1));
pbaspect([2*XBinLim 2*YBinLim 1])
set(gca, 'color', 'none');
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
end