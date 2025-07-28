v = VideoReader('C:\Users\clark\Joe\Github\Gap_Crossing\Data\All_Raw_Videos\2022-06-26 11-51-52_IsoD1_lighting_top_bottom_Exp1.mov');
% vW = VideoWriter('HeadAndCentroidTrackedFly');

% open(vW);

flyCounter = 1;

figure

for flipCounter = 2:4
data = WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(flipCounter);

VerticalTheta = ...
    mod((mod(-1*WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(flipCounter).Orientation(1) + 90,180)-90),180)-90;
            
% MajorAxisToFocus = 2;
headAmbigThresh = 5;
velThresh = 5;

MajorAxisToFocus = 2*data.MajorAxisLength(1)/sqrt(data.MajorAxisLength(1)^2 - data.MinorAxisLength(1)^2);

numFrames = length(data.AbsoluteTime);
avgWindow = 3;
avgVelY = mean(data.VelY(1:avgWindow));
if avgVelY < 0
    data.HeadX(1) = data.CentroidX(1) + data.MajorAxisLength(1)/MajorAxisToFocus*sind(VerticalTheta);
    data.HeadY(1) = data.CentroidY(1) - data.MajorAxisLength(1)/MajorAxisToFocus*cosd(VerticalTheta);
else
    data.HeadX(1) = data.CentroidX(1) - data.MajorAxisLength(1)/MajorAxisToFocus*sind(VerticalTheta);
    data.HeadY(1) = data.CentroidY(1) + data.MajorAxisLength(1)/MajorAxisToFocus*cosd(VerticalTheta);
end

for i = 2:numFrames-1
% for i = 45:50
    VerticalTheta = ...
        mod((mod(-1*WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips(flipCounter).Orientation(i) + 90,180)-90),180)-90;
    frame = read(v,data.AbsoluteTime(i));
    MajorAxisToFocus = 2*data.MajorAxisLength(i)/sqrt(data.MajorAxisLength(i)^2 - data.MinorAxisLength(i)^2);
    tempHeadX1 = data.CentroidX(i) + data.MajorAxisLength(i)/MajorAxisToFocus*sind(VerticalTheta);
    tempHeadY1 = data.CentroidY(i) - data.MajorAxisLength(i)/MajorAxisToFocus*cosd(VerticalTheta);
    tempHeadX2 = data.CentroidX(i) - data.MajorAxisLength(i)/MajorAxisToFocus*sind(VerticalTheta);
    tempHeadY2 = data.CentroidY(i) + data.MajorAxisLength(i)/MajorAxisToFocus*cosd(VerticalTheta);
    predHeadX = data.HeadX(i-1)+data.VelX(i-1);
    predHeadY = data.HeadY(i-1)+data.VelY(i-1);
    if abs(data.VelY(i-1)) > velThresh
        if data.VelY(i-1) < 0
            data.HeadX(i) = tempHeadX1;
            data.HeadY(i) = tempHeadY1;
        else
            data.HeadX(i) = tempHeadX2;
            data.HeadY(i) = tempHeadY2;
        end
    elseif sqrt((tempHeadX1-predHeadX)^2 + (tempHeadY1-predHeadY)^2) - ...
       sqrt((tempHeadX2-predHeadX)^2 + (tempHeadY2-predHeadY)^2) < ...
            headAmbigThresh
        data.HeadX(i) = tempHeadX1;
        data.HeadY(i) = tempHeadY1;
    else
        data.HeadX(i) = tempHeadX2;
        data.HeadY(i) = tempHeadY2;
    end
    % imshow(frame);
    imshow(frame(round(data.CentroidY(i))-50:round(data.CentroidY(i))+50,round(data.CentroidX(i))-50:round(data.CentroidX(i))+50,1))
    hold on
    plot(50,...
         50,...
         'r*');
    plot(data.HeadX(i)-round(data.CentroidX(i))+51,data.HeadY(i)-round(data.CentroidY(i))+51,'b*')
    % plot(data.CentroidX(i),...
    %      data.CentroidY(i),...
    %      'r*');
    % plot(data.HeadX(i),data.HeadY(i),'b*')
    hold off
    pause(1/30)
    % F = getframe;
    % writeVideo(vW,F);
end
end

% close(vW);



% IsoD1MatrixX = IsoD1MatrixX + (IsoD1MatrixMajorAxisLength*MMPerPix/4).*sind(IsoD1MatrixTheta);
% IsoD1MatrixY = IsoD1MatrixY + (IsoD1MatrixMajorAxisLength*MMPerPix/4).*cosd(IsoD1MatrixTheta);