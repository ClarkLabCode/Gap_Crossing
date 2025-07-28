close all

horizShift = 0.15;

% Plot one random good trajectory over time

flyNum = randi(size(IsoD1MatrixX,3)); % 2,19,15
flipNum = randi(size(IsoD1MatrixX,2)); % 18,56,3
flyNum = 2; % Selected this by sifting through random trajectories
flipNum = 18; % Selected this by sifting through random trajectories
oddEven = 1;

f_traj = figure;
f_traj.Position(2) = f_traj.Position(2)/1.75;
f_traj.Position(3:4) = round(1.65*f_traj.Position(3:4));
f_traj.Position(3) = 200;
f_traj.Position(1:2) = [682,302];

XBinLim = 10;

hold on
for i = 1:(size(IsoD1MatrixX,1)-1)
    plot(IsoD1MatrixX(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixY(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
end
xlim([-XBinLim,XBinLim])
ylim([-25,25])
title(['Flip Num ', num2str(flipNum), ', Fly Num ', num2str(flyNum)])
pbaspect([2*XBinLim, 2*25, 1])

% a = plot([0;StraightenedSkelX(end/2+1:end,1);0],...
%          [StraightenedSkelY(19,1);StraightenedSkelY(end/2+1:end,1);StraightenedSkelY(36,1)],'k');
a = plot([StraightenedSkelX(:,oddEven);StraightenedSkelX(1,oddEven)]-horizShift,...
         [StraightenedSkelY(:,oddEven);StraightenedSkelY(1,oddEven)],'k');
uistack(a,'top');
uistack(a,'top');
hold off
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
traj = gca;

f_yt = figure;
f_yt.Position(2:4) = f_traj.Position(2:4);
f_yt.Position(1) = f_traj.Position(1)-f_traj.Position(3);
f_yt.Position(3) = 200;

hold on
for i = 1:(size(IsoD1MatrixX,1)-1)
    plot(IsoD1MatrixTime(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixY(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
end
hold off
xlabel('Time (sec)');
ylabel('Y Position (mm)');
ylim([-25,25])
yt = gca;
pbaspect([20, 2*40, 1])

% f_velyt = figure;
% f_velyt.Position(2:4) = f_traj.Position(2:4);
% f_velyt.Position(1) = f_yt.Position(1)-f_traj.Position(3);
% f_velyt.Position(3) = 200;
% 
% hold on
% for i = 1:(size(IsoD1MatrixX,1)-1)
%     plot(IsoD1MatrixTime(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixVelY(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
% end
% yline(0,'k--');
% hold off
% xlabel('Time (sec)');
% ylabel('Y Velocity (mm/sec)');
% ylim([-40,40])
% xlim([-XBinLim,XBinLim])
% velyt = gca;
% pbaspect([2*XBinLim, 2*40, 1])

f_xt = figure;
f_xt.Position(1,3) = f_traj.Position(1,3);
f_xt.Position(4) = f_traj.Position(4)/3;
f_xt.Position(2) = f_traj.Position(2)-f_xt.Position(4);
f_xt.Position(3) = 120;
f_xt.Position(1) = f_traj.Position(1) + f_traj.Position(3)/2 - f_xt.Position(3)/2;

hold on
for i = 1:(size(IsoD1MatrixX,1)-1)
    plot(IsoD1MatrixX(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixTime(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
end
hold off
xlabel('X Position (mm)');
ylabel('Time (sec)');
xlim([-XBinLim,XBinLim])
xt = gca;
% pbaspect([10, 8, 1])

% f_velxt = figure;
% f_velxt.Position(1,3) = f_traj.Position(1,3);
% f_velxt.Position(4) = f_traj.Position(4)/3;
% f_velxt.Position(2) = f_traj.Position(2)-f_xt.Position(4);
% f_velxt.Position(3) = 120;
% f_velxt.Position(1) = f_traj.Position(1) + f_traj.Position(3)/2 - f_xt.Position(3)/2;
% 
% hold on
% for i = 1:(size(IsoD1MatrixX,1)-1)
%     plot(IsoD1MatrixVelX(i:i+1,flipNum,flyNum,oddEven), IsoD1MatrixTime(i:i+1,flipNum,flyNum,oddEven),'color',[1-(i-1)/(size(IsoD1MatrixX,1)-1) 0 (i-1)/(size(IsoD1MatrixX,1)-1)], 'LineWidth', 2)
% end
% xline(0,'k--');
% hold off
% xlabel('X Velocity (mm/sec)');
% ylabel('Time (sec)');
% ylim([-XBinLim,XBinLim])
% xlim([-25,25])
% velxt = gca;
% 
% 
% linkaxes([traj, yt],'y');
% linkaxes([traj, xt],'x');