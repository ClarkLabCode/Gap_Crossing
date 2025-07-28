% TRACKHEAD Estimates the position of the fly's head for each fly at all times
%  
%  Takes all data in BehavData of a combined WS and adds two new fields for
%  the x and y positions of the head (HeadX and HeadY). This is done by
%  computing which of the two foci of the fitted ellipse are most likely to
%  be the head. The algorithm used first tries to choose the focus in the
%  direction of motion of the centroid if the fly is moving sufficiently
%  fast. If not, the focus chosen is based on its distance relative to the
%  focus chosen in the previous frame. By always giving priority to the
%  direction of motion of the fly, this algorithm minimizes the number of
%  frames in which the wrong focus may be misassigned. In the event that
%  the fly is not moving sufficiently quickly and the distance of both foci
%  are sufficiently close to the previous frame's chosen focus, preference
%  is given to the focus that is higher along the corridor. This choice is
%  made to ensure consistency in upward crossing events which are the only
%  events we generally analyze. Note that the units and coordinate system
%  of HeadX and HeadY generated in this are the same as CentroidX and
%  CentroidY (i.e. - the origin is the top left of the cassette and all
%  units are in pixels and frames).
%  
%  Meant to be used with the output of CombineExpStructs and then used as 
%  the input for AlignCorridors

function WS = TrackHead(WS)

% Set some parameters that are used to optimize the fitting of the head
headAmbigThresh = 5;    % The minimum distance in pix with which the algorithm cannot choose a preferred focus
velThresh = 5;          % The minimum velocity in pix/frame of the centroid for which foci will be chosen by velocity
avgWindow = 3;          % The number of frames used in the initial choice of the first focus

% Loop through every fly
for flyCounter = 1:length(WS.FlipBinnedFlyStruct)
    % Loop through odd and even flips
    for oddEven = 1:2
        % When doing odd flips, load in OddFlips data
        if oddEven == 1
            data = WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips;
        % When doing even flips, load in EvenFlips data
        else
            data = WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips;
        end
        % Loop through each flip over the course of the experiment
        for flipCounter = 1:length(data)
            % If there was no behavior in this flip, skip it
            if length(data(flipCounter).CentroidX) == 1
                continue
            end
            % Compute the scaling factor to convert from major axis length to distance of focus
            MajorAxisToFocus = 2*data(flipCounter).MajorAxisLength(1)/sqrt(data(flipCounter).MajorAxisLength(1)^2 - data(flipCounter).MinorAxisLength(1)^2);
            
            % Compute the angle relative to the vertical rather than the horizontal
            % For a more careful description of why this formulation is
            % required, please see the documentation in AlignCorridors.m
            VerticalTheta = ...
                mod((mod(-1*data(flipCounter).Orientation(1) + 90,180)-90),180)-90;

            % Compute the average y velocity over the first avgWindow frames
            avgVelY = mean(data(flipCounter).VelY(1:avgWindow));

            % If the fly is moving up (which is avgVelY < 0 since the
            % origin is at the top left of the frame), then choose the
            % focus of the ellipse that is higher in the frame
            if avgVelY < 0
                data(flipCounter).HeadX(1) = data(flipCounter).CentroidX(1) + data(flipCounter).MajorAxisLength(1)/MajorAxisToFocus*sind(VerticalTheta);
                data(flipCounter).HeadY(1) = data(flipCounter).CentroidY(1) - data(flipCounter).MajorAxisLength(1)/MajorAxisToFocus*cosd(VerticalTheta);
            % If the fly is moving down, then choose the focus of the
            % ellipse that is lower in the frame
            else
                data(flipCounter).HeadX(1) = data(flipCounter).CentroidX(1) - data(flipCounter).MajorAxisLength(1)/MajorAxisToFocus*sind(VerticalTheta);
                data(flipCounter).HeadY(1) = data(flipCounter).CentroidY(1) + data(flipCounter).MajorAxisLength(1)/MajorAxisToFocus*cosd(VerticalTheta);
            end
            
            % Compute the number of frames in the given flip
            numFrames = length(data(flipCounter).AbsoluteTime);

            % Loop through every frame of the flip but start at the second
            % one since the first frame is initialized above and only go to
            % the second to last frame since we use velocity which is
            % ill-defined in the last frame
            for frameCounter = 2:numFrames-1
                % Compute the angle relative to the vertical rather than the horizontal
                % Again, for a more careful description of why this formulation is
                % required, please see the documentation in AlignCorridors.m
                VerticalTheta = ...
                    mod((mod(-1*data(flipCounter).Orientation(frameCounter) + 90,180)-90),180)-90;

                % Compute the scaling factor to convert from major axis length to distance of focus
                MajorAxisToFocus = 2*data(flipCounter).MajorAxisLength(frameCounter)/sqrt(data(flipCounter).MajorAxisLength(frameCounter)^2 - data(flipCounter).MinorAxisLength(frameCounter)^2);

                % Compute the x and y position of both foci
                % [X1,Y1] corresponds to the higher focus
                tempHeadX1 = data(flipCounter).CentroidX(frameCounter) + data(flipCounter).MajorAxisLength(frameCounter)/MajorAxisToFocus*sind(VerticalTheta);
                tempHeadY1 = data(flipCounter).CentroidY(frameCounter) - data(flipCounter).MajorAxisLength(frameCounter)/MajorAxisToFocus*cosd(VerticalTheta);
                % [X2,Y2] corresponds to the lower focus
                tempHeadX2 = data(flipCounter).CentroidX(frameCounter) - data(flipCounter).MajorAxisLength(frameCounter)/MajorAxisToFocus*sind(VerticalTheta);
                tempHeadY2 = data(flipCounter).CentroidY(frameCounter) + data(flipCounter).MajorAxisLength(frameCounter)/MajorAxisToFocus*cosd(VerticalTheta);
                
                % Compute where the focus is expected to be based on the
                % trajectory of the centroid and the focus location in the
                % previous frame
                predHeadX = data(flipCounter).HeadX(frameCounter-1)+data(flipCounter).VelX(frameCounter-1);
                predHeadY = data(flipCounter).HeadY(frameCounter-1)+data(flipCounter).VelY(frameCounter-1);

                % Check first if the fly is moving sufficiently fast
                if abs(data(flipCounter).VelY(frameCounter-1)) > velThresh
                    % If it is, check which direction the fly is moving
                    % If it's moving upwards, choose the higher focus
                    if data(flipCounter).VelY(frameCounter-1) < 0   % Remember that the origin is at the top left of the frame
                        data(flipCounter).HeadX(frameCounter) = tempHeadX1;
                        data(flipCounter).HeadY(frameCounter) = tempHeadY1;
                    % If it's moving downwards, choose the lower focus
                    else
                        data(flipCounter).HeadX(frameCounter) = tempHeadX2;
                        data(flipCounter).HeadY(frameCounter) = tempHeadY2;
                    end
                % If the fly is not moving sufficiently fast, then we rely
                % on the distance of the two foci from the focus of the
                % previous frame
                % If the higher focus is closer to the predicted focus than
                % the lower one, or if it is sufficiently ambiguous
                % (controlled by headAmbigThresh), then we choose the
                % higher focus
                elseif sqrt((tempHeadX1-predHeadX)^2 + (tempHeadY1-predHeadY)^2) - ...
                   sqrt((tempHeadX2-predHeadX)^2 + (tempHeadY2-predHeadY)^2) < ...
                        headAmbigThresh
                    data(flipCounter).HeadX(frameCounter) = tempHeadX1;
                    data(flipCounter).HeadY(frameCounter) = tempHeadY1;
                % If the lower focus is closer to the predicted focus than
                % the higher one, then we choose the lower focus
                else
                    data(flipCounter).HeadX(frameCounter) = tempHeadX2;
                    data(flipCounter).HeadY(frameCounter) = tempHeadY2;
                end
            end
        end

        % Once we've gone through and done this for all the data, we update
        % WS with the updated fields for both odd and even flips
        % Odd flips
        if oddEven == 1
            WS.FlipBinnedFlyStruct(flyCounter).BehavData.OddFlips = data;
        % Even flips
        else
            WS.FlipBinnedFlyStruct(flyCounter).BehavData.EvenFlips = data;
        end
    end
end

end