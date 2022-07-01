function WS = CorridorSkeletonFinder(WS)

% Port in the relevant fields from WS
NumGaps                         = WS.NumGaps;
NumCorridors                    = WS.NumCorridors;
inputFileName                   = WS.inputFileName;
flipRate                        = WS.flipRate;
ExampleSkeletonCoordsFilePath   = WS.ExampleSkeletonCoordsFilePath;
ExampleSkeletonImgFilePath      = WS.ExampleSkeletonImgFilePath;

% Create video reader object
reader = VideoReader(['..\..\Data\All_Raw_Videos',inputFileName]);

% Initialize necessary matrices
x = zeros(2*NumGaps + 6,2);
y = zeros(2*NumGaps + 6,2);
Skeleton_x = zeros(NumGaps*8 + 4,2,NumCorridors);
Skeleton_y = zeros(NumGaps*8 + 4,2,NumCorridors);

% Do a loop for odd flips (1st flip) and then one for even flips (2nd flip)
for oddEven = 1:2
    % Initialize SkeletonError to have an "error" so that the while loop
    % executes, allowing the user to generate the skeleton
    SkeletonError = 'Yes';
    while strcmp(SkeletonError,'Yes')
        % Create the figure that will have its points labeled by the user
        label_fig = figure();

        % Read in frame 1 for odd flips and frame 1+(flipRate+1)*FrameRate for even flips
        imshow(read(reader,1+(oddEven-1)*(flipRate+1)*reader.FrameRate))

        % Set the location of this figure such that it is compatible with the
        % example figure that guides the user
        label_fig.Units = 'Normalized';
        label_fig.Position = [0.32 0.15 0.66 0.8];

        % Save the location of smallest gap in the odd (1st) flip
        if oddEven == 1
            % Ask user for the location of the smallest gap
            LocOfSmallestGapInFlip = ...
                questdlg('Where is the smallest gap?', 'Smallest Gap', 'Top', 'Bottom', 'Top');
            LocOfSmallestGapInOddFlip = LocOfSmallestGapInFlip;

        % In even (2nd) flip, location of gap is opposite of odd (1st)
        else
            % If it was top in odd flip, it is bottom in even flip
            if strcmp(LocOfSmallestGapInOddFlip, 'Top')
                LocOfSmallestGapInFlip = 'Bottom';
            % If it was bottom in odd flip, it is top in even flip
            else
                LocOfSmallestGapInFlip = 'Top';
            end
        end

        % Load in the coordinates of the example figure
        Example_Coords = load(['..\..\',ExampleSkeletonCoordsFilePath]);

        % Create the figure that already has its points labeled to guide the user
        example_fig = figure();
        frame = imread(['..\..\',ExampleSkeletonImgFilePath]);
        frameSize = size(frame);

        % If the video has the smallest gap at the bottom, load in the example with
        % the smallest gap also at the bottom
        if strcmp(LocOfSmallestGapInFlip,'Bottom')
            imshow(frame((end+1)-(1:end),:,1));
            x_example = Example_Coords.x_example;
            y_example = frameSize(1)-Example_Coords.y_example;
            legend_x = Example_Coords.legend_x;
            legend_y = Example_Coords.legend_y;
        % If the video has the smallest gap at the top, load in the example with
        % the smallest gap also at the top
        else
            imshow(frame(:,:,1));
            x_example = Example_Coords.x_example;
            y_example = Example_Coords.y_example;
            legend_x = Example_Coords.legend_x;
            legend_y = Example_Coords.legend_y;
        end

        % Set the location of the example figure to be adjacent to the one that the
        % user is labeling
        example_fig.Units = 'Normalized';
        example_fig.Position = [0.01 0.15 0.3 0.8];

        % Overlay the markers of the labeled points in the example figure
        figure(example_fig)
        hold on
        plot(x_example,y_example,'ro','MarkerFaceColor','r');
        for pointCounter = 1:length(x_example)
            text(x_example(pointCounter)+5,y_example(pointCounter),num2str(pointCounter),'Color','r','FontSize',14);
        end
        text(legend_x(1),legend_y(1),'Left-most Corridor','HorizontalAlignment',"center",'Color','b','FontSize',14);
        text(legend_x(2),legend_y(1),'Right-most Corridor','HorizontalAlignment',"center",'Color','b','FontSize',14);
        hold off

        % Give the user an opportunity to cancel if they made a mistake on
        % choosing top or bottom for the location of the gap
        GapError = ...
            questdlg(['Are the gaps in the example figure and the figure ', ...
                      'you are labeling in the same orientation?'], ...
                      'Gap orientation verification', 'Yes', 'No', 'No');
        % If the user made a mistake, cancel out and report which file the
        % error occurred on
        if strcmp(GapError,'No')
            close all
            error(['User indicated that the gap orientation was wrong in file: ', inputFileName]) 
        end

        % Have the user input one point at a time
        for pointCounter = 1:length(x_example)
            % Prompt user in command window with instructions
            fprintf(['Label point ', num2str(pointCounter), ...
                ' starting from the left-most corridor and finishing with the', ...
                ' right-most corridor that contains a fly.\n']);

            % Pull up label figure and accept one point at a time
            figure(label_fig)

            % Snippet for allowing to zoom
            zoom on;
            figure(gcf) % Make current figure the active window
            % Give user instructions on how to use the zoom feature
            fprintf('Use the mouse wheel to zoom in and out or select the zoom region by dragging the mouse.\n');
            fprintf('Once you are sufficiently zoomed in, press any keyboard button to enable pixel selection.\n');
            pause() % you can zoom with your mouse and when your image is okay, you press any key
            zoom off; % to escape the zoom mode
            % Save the coordinates selected
            [x_temp,y_temp] = ginput(1);
            zoom out;

            % Round the results to put it into pixels
            x(pointCounter, oddEven) = round(x_temp);
            y(pointCounter, oddEven) = round(y_temp);

            % Mark point so the user sees which ones have been marked
            hold on
            plot(x(pointCounter, oddEven),y(pointCounter, oddEven),'go','MarkerFaceColor','g');
            text(x(pointCounter, oddEven)+5,y(pointCounter, oddEven),num2str(pointCounter),...
                'Color','b','FontSize',14);
            hold off

            % Change the color of the point from red to green in the example image
            % once the user has marked that point
            figure(example_fig)
            hold on
            plot(x_example(pointCounter),y_example(pointCounter),'go','MarkerFaceColor','g');
            text(x_example(pointCounter)+5,y_example(pointCounter),num2str(pointCounter),...
                'Color','b','FontSize',14);
            hold off
        end

        % Now compute the skeleton coordinates based on the points selected
        Skeleton_x([1,2,3,5,6,9,10,13,14,17,18,34,35],oddEven,1) = ...
            x([1,2,11,3,4,5,6,7,8,9,10,13,12],oddEven);
        Skeleton_y([1,2,3,5,6,9,10,13,14,17,18,34,35],oddEven,1) = ...
            y([1,2,11,3,4,5,6,7,8,9,10,13,12],oddEven);
        Skeleton_x([4,7,8,11,12,15,16],oddEven,1) = ...
            x([3,4,5,6,7,8,9],oddEven) - (x(2,oddEven) - x(11,oddEven));
        Skeleton_y([4,7,8,11,12,15,16],oddEven,1) = ...
            y([3,4,5,6,7,8,9],oddEven) - (y(2,oddEven) - y(11,oddEven));
        Skeleton_x(19,oddEven,1) = ...
            x(10,oddEven) + (x(2,oddEven) - x(11,oddEven));
        Skeleton_y(19,oddEven,1) = ...
            y(10,oddEven) + (y(2,oddEven) - y(11,oddEven));
        Skeleton_x([20,23,24,27,28,31,32,35,36],oddEven,1) = ...
            x([9,8,7,6,5,4,3,2,1],oddEven) - (x(2,oddEven) - x(12,oddEven));
        Skeleton_y([20,23,24,27,28,31,32,35,36],oddEven,1) = ...
            y([9,8,7,6,5,4,3,2,1],oddEven) - (y(2,oddEven) - y(12,oddEven));
        Skeleton_x([21,22,25,26,29,30,33,34],oddEven,1) = ...
            Skeleton_x([20,23,24,27,28,31,32,35],oddEven,1) - (x(12,oddEven) - x(13,oddEven));
        Skeleton_y([21,22,25,26,29,30,33,34],oddEven,1) = ...
            Skeleton_y([20,23,24,27,28,31,32,35],oddEven,1) - (y(12,oddEven) - y(13,oddEven));

        % Now overlay the skeletons over the video frames
        figure(label_fig)
        hold on
        plot(Skeleton_x(:,oddEven,1),Skeleton_y(:,oddEven,1),'r*');
        for corrCounter = 1:NumCorridors-1
            figure(label_fig)
            hold on
            Skeleton_x(:,oddEven,corrCounter+1) = ...
                Skeleton_x(:,oddEven,corrCounter) + (1/(NumCorridors-1))*(x(14,oddEven) - x(1,oddEven));
            Skeleton_y(:,oddEven,corrCounter+1) = ...
                Skeleton_y(:,oddEven,corrCounter) + (1/(NumCorridors-1))*(y(14,oddEven) - y(1,oddEven));
            plot(Skeleton_x(:,oddEven,corrCounter+1),Skeleton_y(:,oddEven,corrCounter+1),'r*');
        end

        % Give user opportunity to check and input if an error happened with the skeleton
        % If the user does not indicate an error, this terminates the while
        % loop and proceeds the function
        SkeletonError = questdlg('Is there an error in the skeleton?', 'Skeleton Check', 'Yes', 'No', 'Yes');
            
        % Close the figures
        close(example_fig);
        close(label_fig);
    end
    
end

% Update the fields in WS
WS.LocOfSmallestGapInOddFlip = LocOfSmallestGapInOddFlip;
WS.Skeleton_x = Skeleton_x;
WS.Skeleton_y = Skeleton_y;

end
