[file,path] = uigetfile('*.mat','Choose desired WS file','C:\Users\clark\Joe\Gap_Crossing\Data');
absFilePath = fullfile(path,file);

loadedWS = open(absFilePath);

data = loadedWS.WS.FlipBinnedFlyStruct;

NUM_GAPS = 4;
NUM_FLIES = size(data,2);

for ii = 1:NUM_FLIES
    flyDataOdd = data(ii).AlignedData.OddFlips;
    
    for kk = 1:size(flyDataOdd, 2)
        flip_info = flyDataOdd(kk).UniqCompID;
        flip_events_seq = zeros(1,length(flyDataOdd(kk).CompID));
        
        for jj = 1:NUM_GAPS
            %possible initial events
            upPCross = [2*jj-1, 2*jj, 2*jj+1];
            upGCross = [2*jj-1, 2*jj, 2*jj+1];
            upCirc = [2*jj-1, 2*jj, 2*NUM_GAPS+jj+1];
            upRet = [2*jj-1, 2*jj, 2*jj-1];
            
            downPCross = [2*jj+1, 2*jj, 2*jj-1];
            downGCross = [2*jj+1, 2*jj, 2*jj-1];
            downCirc = [2*jj+1, 2*jj, 2*NUM_GAPS+jj+1];
            downRet = [2*jj+1, 2*jj, 2*jj+1];
            
 
            master_list = {upPCross upGCross upCirc upRet;...
                           downPCross downGCross downCirc downRet};
        
            %EVENT CODE:
                %Hundreds digit: 1-PC, 2-GC, 3-Circ, 4-Ret
                %Tens digit: Gap no.
                %Ones digit: 0-Down, 1-Up
            for mm = 1:size(master_list,1)
                for nn = 1:size(master_list, 2)
                    %look for chamber sequence in UniqCompID to return
                    %first index for occurrences
                    occurrences_indicies = strfind(flip_info, master_list{mm, nn});
                    
                    for pp = 1:length(occurrences_indicies)
                        occurrence_frame = flyDataOdd(kk).UniqCompIDIndex(occurrences_indicies(pp)+1);
                        if nn == 1 && mm == 1
                            if flyDataOdd(kk).RenormUpGlassCrossProb(jj).GapID(pp) < 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        elseif nn == 2 && mm == 1
                            if flyDataOdd(kk).RenormUpGlassCrossProb(jj).GapID(pp) > 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        elseif nn == 1 && mm == 2
                            if flyDataOdd(kk).RenormDownGlassCrossProb(jj).GapID(pp) < 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        elseif nn == 2 && mm == 2
                            if flyDataOdd(kk).RenormDownGlassCrossProb(jj).GapID(pp) > 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        else 
                            flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm;
                        end
                    end
                    
                    
                    %occurrencces_abs_time = flyDataOdd(kk).AbsoluteTime(occurrences_frames);
                end
            end
            
        end
        data(ii).AlignedData.OddFlips(kk).EventSeq = flip_events_seq;
        data(ii).AlignedData.OddFlips(kk).EventSeqCondensed = nonzeros(flip_events_seq)';
        eventSeqCheck = [];
        for qq = 1:NUM_GAPS
            if ~isempty(rmmissing(data(ii).AlignedData.OddFlips(kk).RenormUpGlassCrossProb(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.OddFlips(kk).RenormUpGlassCrossProb(qq).GapID)
                    if data(ii).AlignedData.OddFlips(kk).RenormUpGlassCrossProb(qq).GapID(rr) < 0.5
                        eventSeqCheck = [eventSeqCheck, 100+qq*10+1];
                    else
                        eventSeqCheck = [eventSeqCheck, 200+qq*10+1];
                    end
                    
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.OddFlips(kk).UpCircumventionsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.OddFlips(kk).UpCircumventionsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 300+qq*10+1];
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.OddFlips(kk).UpRetreatsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.OddFlips(kk).UpRetreatsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 400+qq*10+1];
                end
            end
            
            if ~isempty(rmmissing(data(ii).AlignedData.OddFlips(kk).RenormDownGlassCrossProb(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.OddFlips(kk).RenormDownGlassCrossProb(qq).GapID)
                    if data(ii).AlignedData.OddFlips(kk).RenormDownGlassCrossProb(qq).GapID(rr) < 0.5
                        eventSeqCheck = [eventSeqCheck, 100+qq*10+0];
                    else
                        eventSeqCheck = [eventSeqCheck, 200+qq*10+0];
                    end
                    
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.OddFlips(kk).DownCircumventionsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.OddFlips(kk).DownCircumventionsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 300+qq*10+0];
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.OddFlips(kk).DownRetreatsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.OddFlips(kk).DownRetreatsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 400+qq*10+0];
                end
            end
            
            data(ii).AlignedData.OddFlips(kk).EventSeqCheck = eventSeqCheck;
        end
    end
end

%Even
for ii = 1:NUM_FLIES
    flyDataEven = data(ii).AlignedData.EvenFlips;
    
    for kk = 1:size(flyDataEven, 2)
        flip_info = flyDataEven(kk).UniqCompID;
        flip_events_seq = zeros(1,length(flyDataEven(kk).CompID));
        
        for jj = 1:NUM_GAPS
            %possible initial events
            upPCross = [2*jj+1, 2*jj, 2*jj-1];
            upGCross = [2*jj+1, 2*jj, 2*jj-1];
            upCirc = [2*jj+1, 2*jj, 2*NUM_GAPS+jj+1];
            upRet = [2*jj+1, 2*jj, 2*jj+1];
            
            downPCross = [2*jj-1, 2*jj, 2*jj+1];
            downGCross = [2*jj-1, 2*jj, 2*jj+1];
            downCirc = [2*jj-1, 2*jj, 2*NUM_GAPS+jj+1];
            downRet = [2*jj-1, 2*jj, 2*jj-1];
            
 
            master_list = {upPCross upGCross upCirc upRet;...
                           downPCross downGCross downCirc downRet};
        
            
            for mm = 1:size(master_list,1)
                for nn = 1:size(master_list, 2)
                    occurrences_indicies = strfind(flip_info, master_list{mm, nn});
                    
                    for pp = 1:length(occurrences_indicies)
                        occurrence_frame = flyDataEven(kk).UniqCompIDIndex(occurrences_indicies(pp)+1);
                        if nn == 1 && mm == 1
                            if flyDataEven(kk).RenormUpGlassCrossProb(jj).GapID(pp) < 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        elseif nn == 2 && mm == 1
                            if flyDataEven(kk).RenormUpGlassCrossProb(jj).GapID(pp) > 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        elseif nn == 1 && mm == 2
                            if flyDataEven(kk).RenormDownGlassCrossProb(jj).GapID(pp) < 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        elseif nn == 2 && mm == 2
                            if flyDataEven(kk).RenormDownGlassCrossProb(jj).GapID(pp) > 0.5
                                flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm; %event type/gap no./updown
                            end
                        else 
                            flip_events_seq(occurrence_frame) = nn*100 + jj*10 + 2-mm;
                        end
                    end
                    
                    
                   
                end
            end
            
        end
        data(ii).AlignedData.EvenFlips(kk).EventSeq = flip_events_seq;
        data(ii).AlignedData.EvenFlips(kk).EventSeqCondensed = nonzeros(flip_events_seq)';
        eventSeqCheck = [];
        for qq = 1:NUM_GAPS
            if ~isempty(rmmissing(data(ii).AlignedData.EvenFlips(kk).RenormUpGlassCrossProb(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.EvenFlips(kk).RenormUpGlassCrossProb(qq).GapID)
                    if data(ii).AlignedData.EvenFlips(kk).RenormUpGlassCrossProb(qq).GapID(rr) < 0.5
                        eventSeqCheck = [eventSeqCheck, 100+qq*10+1];
                    else
                        eventSeqCheck = [eventSeqCheck, 200+qq*10+1];
                    end
                    
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.EvenFlips(kk).UpCircumventionsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.EvenFlips(kk).UpCircumventionsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 300+qq*10+1];
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.EvenFlips(kk).UpRetreatsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.EvenFlips(kk).UpRetreatsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 400+qq*10+1];
                end
            end
            
            if ~isempty(rmmissing(data(ii).AlignedData.EvenFlips(kk).RenormDownGlassCrossProb(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.EvenFlips(kk).RenormDownGlassCrossProb(qq).GapID)
                    if data(ii).AlignedData.EvenFlips(kk).RenormDownGlassCrossProb(qq).GapID(rr) < 0.5
                        eventSeqCheck = [eventSeqCheck, 100+qq*10+0];
                    else
                        eventSeqCheck = [eventSeqCheck, 200+qq*10+0];
                    end
                    
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.EvenFlips(kk).DownCircumventionsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.EvenFlips(kk).DownCircumventionsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 300+qq*10+0];
                end
            end
            if ~isempty(nonzeros(data(ii).AlignedData.EvenFlips(kk).DownRetreatsIndex(qq).GapID))
                for rr = 1:length(data(ii).AlignedData.EvenFlips(kk).DownRetreatsIndex(qq).GapID)
                    eventSeqCheck = [eventSeqCheck, 400+qq*10+0];
                end
            end
            
            data(ii).AlignedData.EvenFlips(kk).EventSeqCheck = eventSeqCheck;
        end
    end
end
for ii = 1:NUM_FLIES
    flyDataOdd = data(ii).AlignedData.OddFlips;
    flyDataEven = data(ii).AlignedData.EvenFlips;
    length_even_data = size(flyDataEven, 2);
    current_fly_seq = [];
    
    for jj = 1:size(flyDataEven, 2)
        flip_info_odd = flyDataOdd(jj).EventSeqCondensed;
        flip_info_even = flyDataEven(jj).EventSeqCondensed;
        current_fly_seq = [current_fly_seq flip_info_odd 502 flip_info_even 502];
    end
    
    data(ii).CompleteEventSeq = current_fly_seq;
end

%CHECK SEQUENCES

for ii =1:NUM_FLIES
    sum(sort([data(ii).AlignedData.EvenFlips(:).EventSeqCondensed]) ~= sort([data(ii).AlignedData.EvenFlips(:).EventSeqCheck]));
end

NUM_GAPS = 4;
NUM_FLIES = size(data,2);
oddEven = 1;

maxFrameWindow = 30;
numFrameBWSamples = 3;
frameListBefore = -1*[-maxFrameWindow:numFrameBWSamples:-1];
frameListAfter = -1*[numFrameBWSamples:numFrameBWSamples:maxFrameWindow];
FRAMES_PRIOR = frameListBefore; % number of frames to go back before the event occurs
NUM_FRAMES = size(FRAMES_PRIOR,2);
FEATURES = ["AlignedFoldedVelY","AlignedFoldedVelX"];
NUM_FEATURES = size(FEATURES, 2);
GAP_OF_INTEREST = 0;
NUM_FOLDS = 10;
num_predictors = NUM_FRAMES * NUM_FEATURES;
events_list = [] ;
predictors_mat = zeros(0, num_predictors); %dims: num events, predictor

for current_fly = 1:NUM_FLIES
    if oddEven == 1
        flyDataOdd = data(current_fly).AlignedData.OddFlips;
        num_flips = size(flyDataOdd, 2);
        for flip = 1:num_flips
            current_flip_event_seq = flyDataOdd(flip).EventSeq;
    
            %get the last two digits representing gap number and direction
            gap_num_direction = mod(current_flip_event_seq, 100);
    
            %filter for events only occuring at the gap of interest (tens digit)
            if GAP_OF_INTEREST == 0
                event_indices = find(mod(gap_num_direction,10) == 1);
            else
                event_indices = find(gap_num_direction >= GAP_OF_INTEREST*10 & gap_num_direction < GAP_OF_INTEREST*10+10 & mod(gap_num_direction,10) == 1);
            end
            event_indices = event_indices((event_indices > max(abs(FRAMES_PRIOR))) & (event_indices < (length(current_flip_event_seq) - max(abs((FRAMES_PRIOR))))));
            % if max(FRAMES_PRIOR) > 0
            %     event_indices = event_indices(event_indices > max(FRAMES_PRIOR));
            % else
            %     event_indices = event_indices(event_indices < (length(current_flip_event_seq) + min(FRAMES_PRIOR)));
            % end
    
            %add event type to list of outcome variables
            events = current_flip_event_seq(event_indices);
            events_list = [events_list; events'];
    
            for curr_event_index = 1:size(event_indices,2)
                %append predictors 
                flip_predictors = [];
                %iterate by feature type
                for feature = 1:NUM_FEATURES
                    current_feature = FEATURES(feature);
                    feature_vals = flyDataOdd(flip).(current_feature);
                    %iterate by how many frames before event
                    for prior_frame = 1:NUM_FRAMES
                        prev_frame = FRAMES_PRIOR(prior_frame);
                        predictor_val = feature_vals(event_indices(curr_event_index)-prev_frame);
                        flip_predictors = [flip_predictors, predictor_val];
    
                    end
                end
                predictors_mat = [predictors_mat; flip_predictors];
            end
        end
    elseif oddEven == 2
        flyDataEven = data(current_fly).AlignedData.EvenFlips;
        num_flips = size(flyDataEven, 2);
        for flip = 1:num_flips
            current_flip_event_seq = flyDataEven(flip).EventSeq;
    
            %get the last two digits representing gap number and direction
            gap_num_direction = mod(current_flip_event_seq, 100);
    
            %filter for events only occuring at the gap of interest (tens digit)
            %AND for only up events (mod 10 must be 1)
            if GAP_OF_INTEREST == 0
                event_indices = find(mod(gap_num_direction,10) == 1);
            else
                event_indices = find(gap_num_direction >= GAP_OF_INTEREST*10 & gap_num_direction < GAP_OF_INTEREST*10+10 & mod(gap_num_direction,10) == 1);
            end
            event_indices = event_indices((event_indices > max(abs(FRAMES_PRIOR))) & (event_indices < (length(current_flip_event_seq) - max(abs((FRAMES_PRIOR))))));
            % if max(FRAMES_PRIOR) > 0
            %     event_indices = event_indices(event_indices > max(FRAMES_PRIOR));
            % else
            %     event_indices = event_indices(event_indices < (length(current_flip_event_seq) + min(FRAMES_PRIOR)));
            % end
    
            %add event type to list of outcome variables
            events = current_flip_event_seq(event_indices);
            events_list = [events_list; events'];
    
            for curr_event_index = 1:size(event_indices,2)
                %append predictors 
                flip_predictors = [];
                %iterate by feature type
                for feature = 1:NUM_FEATURES
                    current_feature = FEATURES(feature);
                    feature_vals = flyDataEven(flip).(current_feature);
                    %iterate by how many frames before event
                    for prior_frame = 1:NUM_FRAMES
                        prev_frame = FRAMES_PRIOR(prior_frame);
                        predictor_val = feature_vals(event_indices(curr_event_index)-prev_frame);
                        flip_predictors = [flip_predictors, predictor_val];
    
                    end
                end
                predictors_mat = [predictors_mat; flip_predictors];
            end
        end
    end
end
orient_list = mod(events_list,10);
gaps_list = floor(mod(events_list,100)/10); % 1 = 1mm, 2 = 1.5, 3 = 2, 4 = 2.5 
events_list = floor(events_list/100); % 1 = Cross, 2 = Glass Circ, 3 = Gap Circ, 4 = Ret

if oddEven == 1
    if FRAMES_PRIOR == frameListBefore
        savedFileName = strjoin([char(erase(file,'.mat')),...
                         '_OddFlips_BeforeEvent',...
                         '_FrameWindow',num2str(maxFrameWindow),...
                         '_FramesBWSamples',num2str(numFrameBWSamples),...
                         '_Predictors',FEATURES],'');
    elseif FRAMES_PRIOR == frameListAfter
        savedFileName = strjoin([char(erase(file,'.mat')),...
                         '_OddFlips_AfterEvent',...
                         '_FrameWindow',num2str(maxFrameWindow),...
                         '_FramesBWSamples',num2str(numFrameBWSamples),...
                         '_Predictors',FEATURES],'');
    end
elseif oddEven == 2
    if FRAMES_PRIOR == frameListBefore
        savedFileName = strjoin([char(erase(file,'.mat')),...
                         '_EvenFlips_BeforeEvent',...
                         '_FrameWindow',num2str(maxFrameWindow),...
                         '_FramesBWSamples',num2str(numFrameBWSamples),...
                         '_Predictors',FEATURES],'');
    elseif FRAMES_PRIOR == frameListAfter
        savedFileName = strjoin([char(erase(file,'.mat')),...
                         '_EvenFlips_AfterEvent',...
                         '_FrameWindow',num2str(maxFrameWindow),...
                         '_FramesBWSamples',num2str(numFrameBWSamples),...
                         '_Predictors',FEATURES],'');
    end
end

savedFileName = strjoin([savedFileName,'.mat'],'');

save(savedFileName,"gaps_list","predictors_mat","events_list")