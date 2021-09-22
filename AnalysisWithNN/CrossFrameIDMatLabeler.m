NumCorridors = 7;
NumGaps = 4;

% Initialize the video reader object to read in clips
reader = VideoReader('2021-07-07 14-11-22_IsoD1_30C_berberine_chloride_40_min_training.mp4');

% Establish the size of the clips
frameWidth = 100;
frameHeight = 75;

% How many events have already been labeled before? Input 1 if none
eventStartNum = 20;

% Go through the loop for each fly at each width for each orientation
% (even/odd flip) for each event number
for eventCounter = eventStartNum:length(CrossFrameIDMat(1,1,1,:))
    i = 0;
    for gapCounter = 1:NumGaps
        for flyCounter = 1:NumCorridors
            % Check that there is indeed this many events for a given fly
            % at a given width for odd flips
            if CrossFrameIDMat(flyCounter,gapCounter,1,eventCounter) ~= 0
                % Create a writer object to save the clip and open it
                writer = VideoWriter(['Fly_', num2str(flyCounter), '_Gap_', num2str(gapCounter), '_Event_', num2str(eventCounter), '_Odd']);
                open(writer);
                % Read in the frames of the clip where the first frame is
                % 10 frames before the crossing event (if there are 10
                % frames available before the event)
                for frameCounter = max((CrossFrameIDMat(flyCounter,gapCounter,1,eventCounter)-10),1):CrossFrameIDMat(flyCounter,gapCounter,2,eventCounter)
                    frame = read(reader, frameCounter);
                    % Only reads in one channel since not RGB
                    frame = frame(:,:,1);
                    % Only read in the pixels of interest
                    frame = frame((y_cents_1(gapCounter,flyCounter)-round(frameHeight/2)):(y_cents_1(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                 (x_cents_1(gapCounter,flyCounter)-round(frameWidth/2)):(x_cents_1(gapCounter,flyCounter)+round(frameWidth/2)));
                    imshow(frame)
                    
                    % Write frame to video
                    writeVideo(writer, frame);
                    % Slow down the playback of video
                    pause(0.05);
                end
                % Close writer object
                close(writer)
                
                % Request user input to proceed and bring back to command window
                pause;
                close;
                
                % Just a way to track progress when labeling
                i = i+1
                
                % Ask user whether: 
                % (a) Non-glass crossing
                % (b) Glass crossing
                % (c) Leaning non-glass crossing
                % (d) Leaning glass crossing
                % (e) Replay
                % If an input other than (a)-(e) is selected, it defaults to replaying
                % Note this script only allows for one replay per clip
                userInput = input('[a] Non-glass crossing\n[b] Glass crossing\n[c] Leaning non-glass crossing\n[d] Leaning glass crossing\n[e] Replay\n','s');
                if userInput == 'a'
                    % If non-glass crossing, label as 1
                    CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = 1;
                elseif userInput == 'b'
                    % If glass crossing, label as -1
                    CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = -1;
                elseif userInput == 'c'
                    % If leaning non-glass crossing, label as 0.5
                    CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = 0.5;
                elseif userInput == 'd'
                    % If leaning glass crossing, label as -0.5
                    CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = -0.5;
                % Else replay by repeating above loop
                else
                    for frameCounter = max((CrossFrameIDMat(flyCounter,gapCounter,1,eventCounter)-10),1):CrossFrameIDMat(flyCounter,gapCounter,2,eventCounter)
                        frame = read(reader, frameCounter);
                        frame = frame(:,:,1);
                        frame = frame((y_cents_1(gapCounter,flyCounter)-round(frameHeight/2)):(y_cents_1(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                     (x_cents_1(gapCounter,flyCounter)-round(frameWidth/2)):(x_cents_1(gapCounter,flyCounter)+round(frameWidth/2)));
                        imshow(frame)
                        % Play back the clip even slower
                        pause(0.1);
                    end
                
                    % Request user input to proceed and bring back to command window
                    pause;
                    close;

                    % Ask user whether:  
                    % (a) Non-glass crossing
                    % (b) Glass crossing
                    % (c) Leaning non-glass crossing
                    % (d) Leaning glass crossing
                    userInput = input('[a] Non-glass crossing\n[b] Glass crossing\n[c] Leaning non-glass crossing\n[d] Leaning glass crossing\n','s');
                    if userInput == 'a'
                        % If non-glass crossing, label as 1
                        CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = 1;
                    elseif userInput == 'b'
                        % If glass crossing, label as -1
                        CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = -1;
                    elseif userInput == 'c'
                        % If leaning non-glass crossing, label as 0.5
                        CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = 0.5;
                        % If no choice is made after replay, auto labeled as leaning glass
                    else
                        % If leaning glass crossing, label as -0.5
                        CrossFrameIDLabels(flyCounter,gapCounter,eventCounter) = -0.5;
                    end
                end
            end
            
            % Check that there is indeed this many events for a given fly
            % at a given width for even flips
            if CrossFrameIDMat(flyCounter,gapCounter+NumGaps,1,eventCounter) ~= 0
                % Create a writer object to save the clip and open it
                writer = VideoWriter(['Fly_', num2str(flyCounter), '_Gap_', num2str(gapCounter), '_Event_', num2str(eventCounter), '_Even']);
                open(writer);
                % Read in the frames of the clip where the first frame is
                % 10 frames before the crossing event (if there are 10
                % frames available before the event)
                for frameCounter = max((CrossFrameIDMat(flyCounter,gapCounter+NumGaps,1,eventCounter)-10),1):CrossFrameIDMat(flyCounter,gapCounter+NumGaps,2,eventCounter)
                    frame = read(reader, frameCounter);
                    % Only reads in one channel since not RGB
                    frame = frame(:,:,1);
                    % Only read in the pixels of interest
                    frame = frame((y_cents_2(gapCounter,flyCounter)-round(frameHeight/2)):(y_cents_2(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                 (x_cents_2(gapCounter,flyCounter)-round(frameWidth/2)):(x_cents_2(gapCounter,flyCounter)+round(frameWidth/2)));
                    imshow(frame)
                    
                    % Write frame to video
                    writeVideo(writer, frame);
                    % Slow down the playback of video
                    pause(0.05);
                end
                % Close writer object
                close(writer)
                
                % Request user input to proceed and bring back to command window
                pause;
                close;
                
                % Ask user whether: 
                % (a) Non-glass crossing
                % (b) Glass crossing
                % (c) Leaning non-glass crossing
                % (d) Leaning glass crossing
                % (e) Replay
                % If an input other than (a)-(e) is selected, it defaults to replaying
                % Note this script only allows for one replay per clip
                userInput = input('[a] Non-glass crossing\n[b] Glass crossing\n[c] Leaning non-glass crossing\n[d] Leaning glass crossing\n[e] Replay\n','s');
                if userInput == 'a'
                    % If non-glass crossing, label as 1
                    CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = 1;
                elseif userInput == 'b'
                    % If glass crossing, label as -1
                    CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = -1;
                elseif userInput == 'c'
                    % If leaning non-glass crossing, label as 0.5
                    CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = 0.5;
                elseif userInput == 'd'
                    % If leaning glass crossing, label as -0.5
                    CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = -0.5;
                % Else replay by repeating above loop
                else
                    for frameCounter = max((CrossFrameIDMat(flyCounter,gapCounter+NumGaps,1,eventCounter)-10),1):CrossFrameIDMat(flyCounter,gapCounter+NumGaps,2,eventCounter)
                        frame = read(reader, frameCounter);
                        frame = frame(:,:,1);
                        frame = frame((y_cents_2(gapCounter,flyCounter)-round(frameHeight/2)):(y_cents_2(gapCounter,flyCounter)+round(frameHeight/2)), ...
                                     (x_cents_2(gapCounter,flyCounter)-round(frameWidth/2)):(x_cents_2(gapCounter,flyCounter)+round(frameWidth/2)));
                        imshow(frame)
                        % Play back the clip even slower
                        pause(0.1);
                    end
                
                    % Request user input to proceed and bring back to command window
                    pause;
                    close;

                    % Ask user whether:  
                    % (a) Non-glass crossing
                    % (b) Glass crossing
                    % (c) Leaning non-glass crossing
                    % (d) Leaning glass crossing
                    userInput = input('[a] Non-glass crossing\n[b] Glass crossing\n[c] Leaning non-glass crossing\n[d] Leaning glass crossing\n','s');
                    if userInput == 'a'
                        % If non-glass crossing, label as 1
                        CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = 1;
                    elseif userInput == 'b'
                        % If glass crossing, label as -1
                        CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = -1;
                    elseif userInput == 'c'
                        % If leaning non-glass crossing, label as 0.5
                        CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = 0.5;
                        % If no choice is made after replay, auto labeled as leaning glass
                    else
                        % If leaning glass crossing, label as -0.5
                        CrossFrameIDLabels(flyCounter,gapCounter+NumGaps,eventCounter) = -0.5;
                    end
                end
            end
        end
    end
    
    % Advance the eventStartNum to save progress
    eventStartNum = eventStartNum + 1;
    
    % Give the user the option to stop labeling for this session
    pauseScript = input('Do you want to take a break from labeling and save your progress? y/n\n','s');
    if pauseScript == 'y'
        break
    end
end