function [ analysis, allWalkingSpeeds , allTurningSpeeds , allForwardSpeeds, processedData] = AnalyzeFicTracDataBlockAverage(fictracPath, numIgnore, buffer, forceFictrac, redoProjections)
% Processes fictracData.dat of the current data if it exists
% Used in PlotActiveInactiveTimeTraces and MaxResponse
% Params not loaded because this isn't used for plotting
% BA 06/10/2025
% block average the fictrac data like in PlotFicTracAxesBlockAverage
% the forward walking was so noisy in LC15 recordings.

% analysis outputs all trial by trial data even though only
% trialTurningSpeed says trial

% trial turning speed
% trial walking speed
% trial fictrac heading
% trial computed heading by integrating turns

% analysis is a cell containing data for each epoch (ex. analysis{3} = the
% data pertaining to the 3rd epoch in your param file)
% processedData is the raw fictrac data table with ordering: 
%         processedData(i,1) = time;
%         processedData(i,2) = epoch;
%         processedData(i,3) = da_x; % rotation in the x-axis (roll)
%         processedData(i,4) = da_y; % rotation in the y-axis (pitch)
%         processedData(i,5) = da_z; % rotation in the z-axis (yaw)
%         processedData(i,6) = da_t; % time between each frame (seconds)
%         processedData(i,7) = int_x; % integrated x position (mm)
%         processedData(i,8) = int_y; % integrated y position (mm)
%         processedData(i,9) = head; % heading direction
%         processedData(i,10) = instMovDir; 
%         processedData(i,11) = instSpeed;
if ~forceFictrac && isfile([fictracPath '\stimulusData\fictracAnalysis.mat']) && isfile([fictracPath '\stimulusData\allWalkingSpeeds.mat']) && isfile([fictracPath '\stimulusData\processedFicData.mat']) && isfile([fictracPath '\stimulusData\allTurningSpeeds.mat'])
    disp('Previous fictrac analysis found, loading that');
    analysis = load([fictracPath '\stimulusData\fictracAnalysis.mat']).analysis;
    allWalkingSpeeds = load([fictracPath '\stimulusData\allWalkingSpeeds.mat']).allWalkingSpeeds;
    allTurningSpeeds = load([fictracPath '\stimulusData\allTurningSpeeds.mat']).allTurningSpeeds;
    processedData = load([fictracPath '\stimulusData\processedFicData.mat']).processedData;
    allForwardSpeeds = load([fictracPath '\stimulusData\allForwardSpeeds.mat']).allForwardSpeeds;
else

    ficDat = readmatrix([fictracPath '\stimulusData\fictracData.dat']); % read ficDat as a matrix that seperates commas
    params = load([fictracPath '\stimulusdata\chosenparams.mat']).params;
    ficDat(isnan(ficDat)) = 0; % remove all the NaNs
    ficTable = readtable([fictracPath '\stimulusData\fictracData.dat'], 'Delimiter', ';', 'ReadVariableNames',0); % open this row by row with the commas still intact for parsing, using table because it was easier to open it this way
    
    % We read from fictrac every frame and write the results to fictracData.dat 
    % Each time a write is made, it is appended with '|epochNumber|frameTime|'
    % This write can take place in random places throughout the data, so we
    % scan ficTable which has each row of data as a single string. 
    % This tells us what epoch and time that row and all the following rows
    % until the next '|epochNumber|frameTime|' marker takes place. (the buffer
    % usually holds ~2 frames of fictrac data at 140Hz fictrac rate)
    for i = 1:height(ficTable)
        if contains(ficTable.Var1{i}, '|')
            c = split(ficTable.Var1{i},'|');
            ficDat(i,27) = str2double(c{2}); % col 27: epoch number
            ficDat(i,1) = str2double(c{3}); % col 1: read time/frameTime
        end
    end
    
    % get rid of the bad read rows with mumbo jumbo when there were no fictrac reads for a certain frame, this was why
    % some epochs had different numbers of rows
    badInd = find(~ismember(ficDat(:,27), 0:max(ficDat(:,27)))); % not integer epochs
    if ~isempty(badInd)
        disp("Bad rows found, check why");
        disp(badInd);
        ficDat(badInd,:) = []; % replace with zeros instead?
    else
        disp("No bad rows in fictracData, looks good!");
    end
    
    % interpolate frames where fictrac lagged and there was no read
    missedIndices = find(~ficDat(:,2));
    if ~isempty(missedIndices)
        disp(['Fictrac missed ', num2str(length(missedIndices)/length(ficDat)*100,2), '% of the frames']);
        for i = missedIndices'
            if i==length(ficDat) % if the last row of data is a misread
                % repeat the previous row
                ficDat(i,2:26) = ficDat(i-1,2:26);
            elseif i == 1 % if first row is misread
                disp('First row misread, repeating second row');
                ficDat(i,2:26) = ficDat(i+1,2:26);
            else
                % average the surrounding rows
                % if the surrounding rows are zeros too then we got bigger
                % problems
                ficDat(i,2:26) = (ficDat(i+1,2:26)+ficDat(i-1,2:26))/2;
            end
        end
    end

    %% Redoing fictrac camera to animal projections
    if redoProjections
        % old c2a_r was [-2.386918 -0.024729 -1.962130];
%         c2a_r = [-2, -0.4, -2]; % manually annotated transformation
%         c2a_r = [-2.230820, -0.362263, -1.801306]; % steeper x 05282025
%         c2a_r = [-1.912576, -0.460059, -2.022398 ];
%         c2a_r = [0,0,0];
%         c2a_r = [-1.58376526,-0.914076,-1.5825045]; % 100 deg e_z
        % c2a_r = [-1.46198606,-1.02288116,-1.46082231]; % 110 deg e_z
        % c2a_r = [-1.385000000000000,-1.325000000000000,-1.150000000000000];
        c2a_r = [-1.5835,-0.9142,-1.5835] ; % measured angles scipy 'XZY' [0,-120*pi/180,-pi/2]
        c2a_R = omegaToMatrix(c2a_r); % Rodriguez rotation matrix
        ficDat(:,7:9) = ficDat(:,3:5)*c2a_R'; % transpose makes it work
    end
    
    %% block average da_xyz
    groupSize = 10;
    buffer = round(buffer/groupSize); % 30 frames at 14Hz
    % Make sure the vector length is a multiple of groupSize
    N = floor(length(ficDat) / groupSize) * groupSize;
    ficDatTrimmed = ficDat(1:N,:);  % Trim to nearest multiple
    
    % Reshape and compute mean
    da_x_avged = reshape(ficDatTrimmed(:,7), groupSize, []);
    da_x_avged = mean(da_x_avged, 1).'; 

    da_y_avged = reshape(ficDatTrimmed(:,8), groupSize, []);
    da_y_avged = mean(da_y_avged, 1).'; 

    da_z_avged = reshape(ficDatTrimmed(:,9), groupSize, []);
    da_z_avged = mean(da_z_avged, 1).'; 
    % use dT_avged for speed calculations
    dT_avged = reshape(ficDatTrimmed(:,25), groupSize, []);
    dT_avged = mean(dT_avged,1).'/1000000000;
    % for the code to calculate epochTrace, still use the dT_summed
    dT_summed = reshape(ficDatTrimmed(:,25), groupSize, []);
    dT_summed = sum(dT_summed,1).'/1000000000;

    % epochVector = reshape(ficDatTrimmed(:,27), groupSize, []);
    % epochModes = nan(size(epochVector, 2),1);
    % 
    % % Compute mode ignoring zeros
    % for i = 1:size(epochVector, 2)
    %     group = epochVector(:, i);
    %     group_nonzero = group(group ~= 0);
    %     if ~isempty(group_nonzero)
    %         epochModes(i) = mode(group_nonzero);
    %     end
    % end

    % find epoch
    timeSampled = cumsum(dT_summed); % changing this to dT_avged instead of dT_summed broke this algo
    epochTrace = nan(length(timeSampled),1);
    for i = 1:length(timeSampled)
        [~, idx] = min(abs(ficDatTrimmed(:,1) - timeSampled(i)));
        epochTrace(i) = ficDatTrimmed(idx,27);
    end

    %% Low Pass filter the rotation values in raw ficDat table (try doing filter and then average)
    filt_order = 6;
    fc = 5; % cutoff frequency
    fs = 14; % fictrac is sampled at 140Hz
    [b,a] = butter(filt_order, fc/(fs/2),'low');
    da_x_avged = filtfilt(b,a,da_x_avged);
    da_y_avged = filtfilt(b,a,da_y_avged);
    da_z_avged = filtfilt(b,a,da_z_avged);

    %% Write to the 14Hz processedData matrix
    sphereRadius = 4.5; % mm
    % initialize processed data array
    % frame time, epoch, dax, day, daz, dt, intx, inty, intHead, instMovDir, instSpeed
    % processedData = zeros(length(epochInds), 12);
    processedData = zeros(length(epochTrace),12);
    
    % fictracData already has the correct frame timings from Q.timing.flipt-Q.timing.t0;
    % processedData(:,1) = ficDat(epochInds,1); % 0, 0.0168, 0.0333,... TODO: closedloop data has 2/60 frame rate?
    processedData(:,1) = timeSampled; % nonzero start
    processedData(:,2) = epochTrace;
    processedData(:,3) = da_x_avged;
    processedData(:,4) = da_y_avged;
    processedData(:,5) = da_z_avged;
    processedData(:,6) = dT_avged;
    processedData(:,11) = sqrt(da_x_avged.^2+da_y_avged.^2+da_z_avged.^2)*sphereRadius./dT_avged;

    % filter out outliers
    outlierInds = processedData(:,11) > 15; % any inst speeds greater than 15mm/s
    processedData(outlierInds,3) = 0; processedData(outlierInds,4) = 0; processedData(outlierInds,5) = 0; processedData(outlierInds,11) = 0;


    %% Output all walking and turning speeds
%     allWalkingSpeeds = filtfilt(b,a,processedData(:,11));
%     allTurningSpeeds = filtfilt(b,a,processedData(:,5)./processedData(:,6)*180/pi);
    % allForwardSpeeds = filtfilt(b,a,processedData(:,4)./processedData(:,6)*sphereRadius);
    allWalkingSpeeds = processedData(:,11);
    allTurningSpeeds = processedData(:,5)./processedData(:,6)*180/pi;
    allForwardSpeeds = processedData(:,4)./processedData(:,6)*sphereRadius;
    % allTurningSpeeds = ficDat(:,9)./ficDat(:,25)*1000000000*180/pi;
    % allForwardSpeeds = ficDat(:,8)./ficDat(:,25)*1000000000*180/pi;
    
    %% calculate turning and walking speeds
    % TODO: add more processedData fields into the analysis cell
    % dv_theta = processedData(:,5)./processedData(:,6)*180/pi;
    
    numEpochs = max(processedData(:,2));
    analysis = cell(numEpochs, 1);
    for i = numIgnore+1:numEpochs % which stim epoch 
        thisEpochInds = find(ismember(processedData(:,2),i)); % find all the row indices corresponding to this epoch
        epEnd = [0; find(diff(thisEpochInds)>25)]; % find the indices of the above matrix that correspond to epochs ending assuming that the interleave is longer than 25 frames
%         if length(epEnd) == 1
%             epEnd = [0];
%         end
        analysis{i}.trialTurningSpeed = cell(length(epEnd), 1); % length(epEnd) = numRepeats of the epoch
        analysis{i}.trialWalkingSpeed = cell(length(epEnd), 1);
        analysis{i}.fictracHeading = cell(length(epEnd),1); % the heading calculated by fictrac
        analysis{i}.compHeading = cell(length(epEnd),1); % the cumulative da_z turning, should match fictracHeading
        % epLength = round(length(thisEpochInds)/length(epEnd)); % assume all trials have the same number of frames...
        epLength = round(params(i).duration/60*14);
        for j = 1:length(epEnd) % looking at each trial
            analysis{i,1}.trialTimeX = linspace(-buffer/14,params(i).duration/60+buffer/14,epLength+2*buffer).'; % -30frames/14Hz to epDuration (s) + 30frames/14Hz
            if j == length(epEnd) % last repeat/trial
                try
                    epda_z = -processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,5);
                    epda_y = -processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,4);
                    epdt = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,6);
                    analysis{i,1}.trialTurningSpeed{j,1} = epda_z./epdt*180/pi; % (deg/s)
                    analysis{i,1}.trialWalkingSpeed{j,1} = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,11); % (mm/s)
                    analysis{i,1}.trialForwardSpeed{j,1} = epda_y./epdt*sphereRadius; % (mm/s)
                    analysis{i,1}.averageWalkingSpeed{j,1} = mean(analysis{i,1}.trialWalkingSpeed{j,1});
                    analysis{i,1}.fictracHeading{j,1} = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,9); % (rads, modulus)
                    analysis{i,1}.compHeading{j,1} = cumsum(epda_z); % (rads, not modulus)
                catch % last stimuli presentation so the interleave isn't played again afterwards so there is no buffer to plot. Concatenate with zeros
                    epda_z = -processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(end),5);
                    epda_y = -processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(end),4);
                    epdt = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(end),6);
                    missing = epLength+2*buffer-length(epda_z); % epLength includes buffer
                    % missing = length(epda_z) - epLength;
                    if missing < 0
                        disp(num2str(missing));
                        warning('calculating last trial of exp wrong');
                    end
                    analysis{i,1}.trialTurningSpeed{j,1} = [epda_z;zeros(missing,1)]./[epdt;ones(missing,1)]*180/pi; % (deg/s)
                    analysis{i,1}.trialWalkingSpeed{j,1} = [processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(end),11);zeros(missing,1)]; % (mm/s)
                    analysis{i,1}.trialForwardSpeed{j,1} = [epda_y;zeros(missing,1)]./[epdt;ones(missing,1)]*sphereRadius; % (mm/s)
                    analysis{i,1}.averageWalkingSpeed{j,1} = mean(analysis{i,1}.trialWalkingSpeed{j,1});
                    analysis{i,1}.fictracHeading{j,1} = [processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(end),9);zeros(missing,1)]; % (rads, modulus)
                    analysis{i,1}.compHeading{j,1} = cumsum([epda_z;zeros(missing,1)]); % (rads, not modulus)
                end
            else
                epda_z = -processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,5);
                epda_y = -processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,4);
                epdt = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,6);
                analysis{i,1}.trialTurningSpeed{j,1} = epda_z./epdt*180/pi;
                analysis{i,1}.trialWalkingSpeed{j,1} = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,11);
                analysis{i,1}.trialForwardSpeed{j,1} = epda_y./epdt*180/pi;
                analysis{i,1}.averageWalkingSpeed{j,1} = mean(analysis{i,1}.trialWalkingSpeed{j,1});
                analysis{i,1}.fictracHeading{j,1} = processedData(thisEpochInds(epEnd(j)+1)-buffer:thisEpochInds(epEnd(j)+1)+epLength+buffer-1,9); % (rads, modulus)
                analysis{i,1}.compHeading{j,1} = cumsum(epda_z); % (rads, not modulus)
                % I previously used epLength = length(epda_z) here assuming
                % there would always be multiple trials. But when there is only
                % one trial, I have to define it above
    %             epLength = length(epda_z); % jank way of doing this but honestly not a bad way to make sure everything is working lol
            end
        end
    end
    %% saving all analysis variables
%     [fictracPath '\stimulusData\fictracAnalysis.mat']) && isfile([fictracPath '\stimulusData\allWalkingSpeeds.mat']) && isfile([fictracPath '\stimulusData\processedFicData.mat'])
    save([fictracPath '\stimulusData\fictracAnalysis.mat'], 'analysis');
    save([fictracPath '\stimulusData\allWalkingSpeeds.mat'], 'allWalkingSpeeds');
    save([fictracPath '\stimulusData\processedFicData.mat'], 'processedData');
    save([fictracPath '\stimulusData\allTurningSpeeds.mat'], 'allTurningSpeeds');
    save([fictracPath '\stimulusData\allForwardSpeeds.mat'], 'allForwardSpeeds')
end

end

