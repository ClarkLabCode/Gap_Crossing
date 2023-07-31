% CORRIDORSKELETONFINDER Prompts user to outline skeleton of cassette via GUI
%
%  When this function is called, the user will first be presented a single
%  frame from the first odd flip and guided through labeling the points
%  necessary to draw out the skeleton of each cassette. Upon completion,
%  the user will be asked to verify their work. If there is a mistake, this
%  step repeats. If no mistake was made, the the user will be presented a
%  single frame from the first even flip and once again guided through
%  labeling. There will be one last check for a mistake. Once everything is
%  confirmed, the skeleton points are saved to WS accordingly.

function WS = CorridorSkeletonFinder(WS)

% Port in the relevant fields from WS
NumGaps                         = WS.NumGaps;
NumCorridors                    = WS.NumCorridors;
inputFileName                   = WS.inputFileName;
flipRate                        = WS.flipRate;
ExampleSkeletonCoordsFilePath   = WS.ExampleSkeletonCoordsFilePath;
ExampleSkeletonImgFilePath      = WS.ExampleSkeletonImgFilePath;

% Create video reader object
reader = VideoReader(['..\0_Raw_Videos\',inputFileName]);

% Initialize necessary matrices
x = zeros(2*NumGaps + 6,2);
y = zeros(2*NumGaps + 6,2);
Skeleton_x = zeros(NumGaps*8 + 4,2,NumCorridors);
Skeleton_y = zeros(NumGaps*8 + 4,2,NumCorridors);

% Establish some values and vectors that will be used to index through the 
% above matrices when performing computations to extract the skeleton
L_Skel = 8*NumGaps + 4; % Number of points in a given orientation's skeleton for a single corridor
L1 = (2*(NumGaps+1) + 3); % Number of skeleton points from operation 1

% Vec1_Skel is the vector of indices in Skeleton_x/y being manipulated in operation 1
Vec1_Skel = zeros(L1,1);
Vec1_Skel(1:4) = [1,2,3,5];
Vec1_Skel((end-2):end) = [L_Skel/2, L_Skel-2, L_Skel-1];
for i = 1:(2*(NumGaps-1))
    Vec1_Skel(i+4) = i + 5 + 2*floor(i/2);
end

% Vec1_x is the vector of indices in x/y being used for manipulations in operation 1
Vec1_x = 1:L1;
Vec1_x = circshift(Vec1_x,1);
Vec1_x(3) = Vec1_x(end-1);
Vec1_x(end-1) = L1;
Vec1_x(1:2) = 1:2;
Vec1_x = Vec1_x';

L2 = 2*NumGaps-1; % Number of skeleton points from operation 2

% Vec2_Skel is the vector of indices in Skeleton_x/y being manipulated in operation 2
Vec2_Skel = zeros(L2,1);
for i = 1:L2
    Vec2_Skel(i) = i + 3 + 2*floor(i/2);
end

% Vec2_x is the vector of indices in x/y being used for manipulations in operation 2 
Vec2_x = 2+(1:L2);
Vec2_x = Vec2_x';

L3 = 2*NumGaps+1; % Number of skeleton points from operation 3

% Vec3_Skel is the vector of indices in Skeleton_x/y being manipulated in operation 3
Vec3_Skel = zeros(L3,1);
Vec3_Skel(end) = L_Skel;
for i = 1:(L3-1)
    Vec3_Skel(i) = i + (4*NumGaps + 2) + 2*floor((i-1)/2);
end

% Vec3_x is the vector of indices in x/y being used for manipulations in operation 3
Vec3_x = 1 + (1:L3);
Vec3_x(1) = Vec3_x(1)-1;
Vec3_x = flip(Vec3_x);
Vec3_x = Vec3_x';

L4 = 2*NumGaps - 1; % Number of skeleton points from operation 4

% Vec4_Skel_a is the vector of indices in Skeleton_x/y being manipulated in operation 4
% Vec4_Skel_b is the vector of indices in Skeleton_x/y being used for manipulations in operation 4
Vec4_Skel_a = zeros(L4,1);
Vec4_Skel_b = zeros(L4,1);
for i = 1:L4
    Vec4_Skel_a(i) = i + (4*NumGaps + 2) + 2*floor((i+1)/2);
    Vec4_Skel_b(i) = i + (4*NumGaps + 3) + 2*floor(i/2);
end

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
        Example_Coords = load(['..\..\..\..\',ExampleSkeletonCoordsFilePath]);

        % Create the figure that already has its points labeled to guide the user
        example_fig = figure();
        frame = imread(['..\..\..\..\',ExampleSkeletonImgFilePath]);
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
        for pointCounter = 1:min(14,size(x,1)) % min 14 is because overlay is for standard 4 gaps
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
        for pointCounter = 1:size(x,1)
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
            if pointCounter <= 14 % Cut off the green points at 14 since example overlay is 4 gaps
                figure(example_fig)
                hold on
                plot(x_example(pointCounter),y_example(pointCounter),'go','MarkerFaceColor','g');
                text(x_example(pointCounter)+5,y_example(pointCounter),num2str(pointCounter),...
                    'Color','b','FontSize',14);
                hold off
            end
        end
        
        % Now compute the skeleton coordinates based on the points selected
        % All points user already selected on skeleton
        Skeleton_x(Vec1_Skel,oddEven,1) = ...
            x(Vec1_x,oddEven);
        Skeleton_y(Vec1_Skel,oddEven,1) = ...
            y(Vec1_x,oddEven);
        % Left side outer gap points not selected by user
        Skeleton_x(Vec2_Skel,oddEven,1) = ...
            x(Vec2_x,oddEven) - (x(2,oddEven) - x((2*(NumGaps+1) + 1),oddEven));
        Skeleton_y(Vec2_Skel,oddEven,1) = ...
            y(Vec2_x,oddEven) - (y(2,oddEven) - y((2*(NumGaps+1) + 1),oddEven));
%         % Top right corner of corridor
%         Skeleton_x((L_Skel/2 + 1),oddEven,1) = ...
%             x(2*(NumGaps+1),oddEven) + (x(2,oddEven) - x((2*(NumGaps+1) + 1),oddEven));
%         Skeleton_y((L_Skel/2 + 1),oddEven,1) = ...
%             y(2*(NumGaps+1),oddEven) + (y(2,oddEven) - y((2*(NumGaps+1) + 1),oddEven));
        % Right side inner gap points not selected by user
        Skeleton_x(Vec3_Skel,oddEven,1) = ...
            x(Vec3_x,oddEven) - (x(2,oddEven) - x((2*(NumGaps+1) + 2),oddEven));
        Skeleton_y(Vec3_Skel,oddEven,1) = ...
            y(Vec3_x,oddEven) - (y(2,oddEven) - y((2*(NumGaps+1) + 2),oddEven));
        % Right side outer gap points not selected by user
        Skeleton_x(Vec4_Skel_a,oddEven,1) = ...
            Skeleton_x(Vec4_Skel_b,oddEven,1) - (x((2*(NumGaps+1) + 2),oddEven) - x((2*(NumGaps+1) + 3),oddEven));
        Skeleton_y(Vec4_Skel_a,oddEven,1) = ...
            Skeleton_y(Vec4_Skel_b,oddEven,1) - (y((2*(NumGaps+1) + 2),oddEven) - y((2*(NumGaps+1) + 3),oddEven));

        % Now overlay the skeletons over the video frames
        figure(label_fig)
        hold on
        plot(Skeleton_x(:,oddEven,1),Skeleton_y(:,oddEven,1),'r*');
        for corrCounter = 1:NumCorridors-1
            figure(label_fig)
            hold on
            Skeleton_x(:,oddEven,corrCounter+1) = ...
                Skeleton_x(:,oddEven,corrCounter) + (1/(NumCorridors-1))*(x((2*(NumGaps+1) + 4),oddEven) - x(1,oddEven));
            Skeleton_y(:,oddEven,corrCounter+1) = ...
                Skeleton_y(:,oddEven,corrCounter) + (1/(NumCorridors-1))*(y((2*(NumGaps+1) + 4),oddEven) - y(1,oddEven));
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
