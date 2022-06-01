%%
% short timescale kinematics of ZF

clear; clc;
%set(0,'DefaultFigureVisible','off');


addpath(genpath('addPathToYourFunctions')); %path to where useful the functions are stored
cd 'addPathToWhereYourScriptsAreLocated'; %path to where the script is stored

%%  set recording settings
datarate_Hz=700;
camscale_px_per_mm=17;
TL=2;% tail length in mm

%%

% load dir
myDir = uigetdir('addPathToWhereYourDataIsStored','Go to your directory!'); %path to where your data is stored
if ~isdir(myDir)
    uiwait(warndlg('The selected directory does not exist'));
    return;
end

filePattern = fullfile(myDir,'*.txt');
myFiles = dir(filePattern); %all txt files

tic; % start timer

for ff = 1:length(myFiles) % go through each file
    
    dirName = myDir;
    fileName = myFiles(ff).name;
    inputName = fullfile(myDir,fileName);

% extract parameters from filename
tmp_str = strsplit(fileName, '_');

% save parameters in strings - FORMAT: Date_Time_ExperimentName_Animal#_Remark_Trial#_.bin
% acquis_date = tmp_str{1, 1}; acquis_time = tmp_str{1, 2}; exp_type    = tmp_str{1, 3}; fish_num    = tmp_str{1, 4}; trial_num   = tmp_str{1, 6};
fish_species=tmp_str{1, 1};dpf=tmp_str{1, 2};animal_number=tmp_str{1, 3};

% acquisition parameters
window_size         = 0; %useful when writing imgs to disk
num_data_categories = 33;

% Read and reorganize the bin file

tmp_data = txt2mat(inputName);

% h        = fopen(inputName);
% tmp_data = fread(h, inf, 'char');
% fclose(h);

%tmp_data = tmp_data(1:(end-mod(size(tmp_data,1), num_data_categories)), 1); % cuts when parts of the data categories are missing at the end
%tmp_data = reshape(tmp_data, [num_data_categories, size(tmp_data, 1)/num_data_categories])';

%% data key
% 1: frame #
% 2: xpos
% 3: ypos
% 6: head yaw angle
% 8: fishConfidence
% 9 and 10 - NA
% 11-20: px value of tail segs
% 21-31: tail seg angles
% 32: frame # stamp
% 33: frames lagging

%% clean bad tracking and replace with nans

% based on fish tracking confidence
fishConfidence=tmp_data(:,8);
histogram(fishConfidence); 
disp("fishConfidence histogram"); %uiwait();

idx_TE1=find(fishConfidence<100);
%tmp_data(idx_TE1,2:end)=NaN; % can be misleading as the erraneous ID is
%often a bright spot; go to tail confidence

% based on tail confidence; sometimes a spot is accidentally detected in the absence of a fish but
% the tail confidence will still remain zero
tailConfidence=sum(tmp_data(:,11:21),2);
histogram(tailConfidence); 
disp("tailConfidence histogram"); %uiwait();

idx_TE2=find(tailConfidence<100);
tmp_data(idx_TE2,2:end)=NaN;

%idx_TE2=find(tailConfidence==0);
%tmp_data(idx_TE2,2:end)=NaN;

tmp_frame=tmp_data(:,1);
frameNum=tmp_frame-tmp_frame(1);
frameDiff=diff(frameNum);
plot(frameDiff); 
disp("freme diff - look at lost frames"); %uiwait();

%% CHECKPOINT
% CHECK FOR LOST FRAMES

% FRAME COUNTER:
% index difference between frames, based on the cameras 24bit frame counter

frame_diff = [0; diff(tmp_data(:, 1))]; 

% CHECK for missing frames
idx_frame  = frameDiff > 1;                         % index of missing frames
idx_lost   = find(idx_frame == 1);                   % first frame in the block of missed frames


% calculate total duration of video

duration = ((frameNum(end,1)*1.4286)/1000)/60;

% prints the above calculated values

fprintf('\nfish species: %s', fish_species);
fprintf('\nage: %s', dpf);
fprintf('\nvideo duration: %2.2f min', duration);
fprintf('\n\nfirst frame in the block of missed frames : number of frames lost\n');
fprintf('\n %d: %d',  [idx_lost, frameDiff(idx_frame)-1].');
fprintf('\nCurrent file is %s\n',fileName);
fprintf('Now reading %.2f',ff); fprintf(' of %.2f files\n',length(myFiles));


% % INSERT nans for lost frames...
% 
% % define anonymous function that inserts (nxm) blocks into (oxm) matrices
insert_blocks = @(insert_block, matrix, n) cat(1,  matrix(1:n-1,:), insert_block, matrix(n:end,:) );

data_raw = tmp_data;

for ii = nnz(idx_frame):-1:1 % starts from the last row in the matrix to keep the indizes for block-insertion
 
    nan_block       = nan(frame_diff(idx_lost(ii)) - 1, num_data_categories);
    nan_block(:, 1) = tmp_data(idx_lost(ii)-1, 1)+1: tmp_data(idx_lost(ii)-1, 1) + frame_diff(idx_lost(ii))-1; % fill the first column of the Nan blocks with frame numbers that were missing
    
    tmp_data        = insert_blocks(nan_block, tmp_data, idx_lost(ii));
    
end

tmp_data(:,1) = tmp_data(:,1) - tmp_data(1,1) + 1; % framecounter starts at 1


%% get basic fish params; plot and check

% xPos and yPos
xPos=tmp_data(:,2);yPos=tmp_data(:,3);
plot(xPos,yPos);

% head yaws
headyaw=tmp_data(:,6);

headyaw_diff = [nan; diff(headyaw)];
headyaw_diff_tmp=headyaw_diff;

% correction of discontinuties when fish turns from 0 to 2*pi or vice versa

for kk = 1: length(headyaw_diff)
    
    if headyaw_diff(kk) > pi % right turn
        headyaw_diff(kk) =  headyaw_diff(kk) - 2*pi;
        
    elseif headyaw_diff(kk) < -pi % left turn
        headyaw_diff(kk) =  2*pi + headyaw_diff(kk);
    end
    
end


% tail segs
tailSeg=tmp_data(1:size(tmp_data,1),11:21);

% tail angles
tailSegAngles=tmp_data(1:size(tmp_data,1),22:31);


%% polynomial filter

xPos_filt = sgolayfilt(xPos,1,21);
yPos_filt = sgolayfilt(yPos,1,21);

plot(xPos_filt,yPos_filt); 
disp("filtered x-y plot"); %uiwait();

plot(xPos);hold on; plot(xPos_filt); hold off; 
disp("raw vs. filtered x-pos"); %uiwait();


%%
dx   = [0; diff(xPos_filt)]; % distance between two consecutive x-coordinates
dy   = [0; diff(yPos_filt)]; % distance between two consecutive y-coordinates


tmp_dist_unfilt           = sqrt(dx.^2 + dy.^2);
tmp_dist_unfilt           = tmp_dist_unfilt./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate_Hz;         % convert to velocity in mm/s


%% temp - to check how completely unfiltered traces look like when velo is cal

dxZ   = [0; diff(xPos)]; % distance between two consecutive x-coordinates
dyZ   = [0; diff(yPos)]; % distance between two consecutive y-coordinates


tmp_dist_unfiltZ           = sqrt(dxZ.^2 + dyZ.^2);
tmp_dist_unfiltZ           = tmp_dist_unfiltZ./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_unfiltRaw            = tmp_dist_unfiltZ.*datarate_Hz;         % convert to velocity in mm/s


%% filtering distances
%change here for DT - in ZF we use various event detection filters wheareas in DT, 
% only a polynomial filter is used to smothen the trace as we don't want to assume anything 
% about the structure of our signal

idx_nan     = isnan(dx);
idx_nan = find(idx_nan==1);

dx(idx_nan) = 0; % for filtering nan values need to be removed 
dy(idx_nan) = 0;

% filters used in the analysis
filterB = ones(1,100)./100; %for event detection 100 ms
filterV = ones(1,30)./30; %precise onset
filterF = ones(1,10)./10; %for escapes 10 ms


%% Event detection
dxB        = filtfilt(filterB, 1, dx);
dyB        = filtfilt(filterB, 1, dy);

tmp_dist_fB = sqrt(dxB.^2 + dyB.^2);     % distance moved between iterations, in pixels
tmp_dist_fB = tmp_dist_fB./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fB  = tmp_dist_fB.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fB(idx_nan) = nan; % re-insert the nan values


%% Onset detection - also follows the signal more closely
dxV        = filtfilt(filterV, 1, dx);
dyV        = filtfilt(filterV, 1, dy);

tmp_dist_fV = sqrt(dxV.^2 + dyV.^2);     % distance moved between iterations, in pixels
tmp_dist_fV = tmp_dist_fV./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fV  = tmp_dist_fV.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fV(idx_nan) = nan; % re-insert the nan values

% re-insert nans in dx dy
dx(idx_nan) = NaN; % for filtering nan values need to be removed 
dy(idx_nan) = NaN;

plot(tmp_vel_unfilt); hold on; plot(tmp_vel_fB); plot(tmp_vel_fV); plot(tmp_vel_unfiltRaw); hold off;
disp("unfilt and broad/narrow filtered traces"); %uiwait();


%% tail angle processing from Mike Orger

tailSegAngles(tailSegAngles<-10)=0;
cumsumtailangles=cumsum(tailSegAngles')';

%automated selection for properly tracked segments; the others are
%discarded in the analysis
pxVal_tail=nanmedian(tailSeg);
numsegs=8;%length(pxVal_tail(pxVal_tail>30))

smoothedCumsumFixedSegmentAngles=cumsumtailangles;

%for the bout detection we smooth the tail curvature to eliminate kinks due
%to tracking noise
for n=2:size(cumsumtailangles,2)-1
    smoothedCumsumFixedSegmentAngles(:,n)=mean(cumsumtailangles(:,n-1:n+1)');
end

%we consider the difference in segment angles, because we want to detect
%tail movement
fixedSegmentAngles = [zeros(size(cumsumtailangles,2),1)'; diff(smoothedCumsumFixedSegmentAngles)];

%just take the number of well tracked segments
realSegmentAngles = fixedSegmentAngles(:,1:numsegs);
disp(1)

%Filter angles
%this parameters set the timescale at which to smooth the tail movement
bcsize=10;
filteredSegmentAngle = realSegmentAngles*0;
for n = 1 : (size(realSegmentAngles,2))
    %smooth the tail movements at the shorter timescale
filteredSegmentAngle(:,n) = conv(realSegmentAngles(:,n),ones(bcsize,1)'/bcsize,'same');
end
disp(2)

%we sum the angle differences down the length of the tail to give
%prominence to regions of continuous curvature in one direction
cumFilteredSegmentAngle = cumsum(filteredSegmentAngle')';

%we sum up the absolute value of this accumulated curvature so bends in
%both directions are considered
superCumSegAngle = cumsum(abs(cumFilteredSegmentAngle)')';

%Calculate tail curvature
%it convolves all the segments in one
%this measure is filtered with a boxcar filter at the tail-beat timescale
bcFilt=10;
tailCurveMeasure = conv(superCumSegAngle(:,end),boxcarf(bcFilt),'same');

%min max filter to smooth data
%min filter removes baseline fluctuations (needs to be well above the bout
%length)
%max filter flattens out beat variations (timescale needs to be adjusted
%according to be well below the interbout length)
smootherTailCurveMeasure = maxFilterFast(tailCurveMeasure(:,end),20) - minFilterFast(tailCurveMeasure(:,end),400);

%fix size bug of smootherTailCurveMeasure
if length(smootherTailCurveMeasure) ~= length(cumsumtailangles)
    smootherTailCurveMeasure(end) = [];
end


%%
%manually selected threshold
allbout=smootherTailCurveMeasure>0.25;
allstart=max([0 diff(allbout)],0);
allend=max([0 -diff(allbout)],0);

%% ID regions to calculate bout associated kinematics

starts=find(allstart);
ends=find(allend);

if starts(1)>ends(1)
ends(1)=[]; % to remove already init swim
end

if starts(end)>ends(end)
starts(end)=[]; % to remove half swim in the end
end


%%
% filter traces

rawTail=cumsumtailangles(:,numsegs);

forTBF=sgolayfilt(cumsumtailangles(:,numsegs-1),3,21); %7th segment
forMaxAngle=abs(sgolayfilt(cumsumtailangles(:,numsegs),2,11));
headyawFilt=sgolayfilt(headyaw_diff,2,11);
plot(cumsumtailangles(:,numsegs)); hold on; plot(forTBF);

% interp small gaps ONLY - 5frames/7.1ms
forTBF_intrp=(interp1gap(forTBF,5,'linear'))';
tmp_dist_fV_intrp=(interp1gap(tmp_dist_fV,5,'linear'))';
tmp_vel_fB_intrp=(interp1gap(tmp_vel_fB,5,'linear'))';
tmp_vel_fV_intrp=(interp1gap(tmp_vel_fV,5,'linear'))';
forMaxAngle_intrp=(interp1gap(forMaxAngle,5,'linear'))';
headyawFilt_intrp=(interp1gap(headyawFilt,5,'linear'))';

% cal accn
forAccn = [0; diff(tmp_vel_fB_intrp)];
forAccnUnfilt = forAccn.*datarate_Hz;
forAccnFilt =sgolayfilt(forAccnUnfilt,2,15);

plot(rawTail); hold on; plot(forTBF_intrp); disp("verify raw 7th segment with filtered one");
%uiwait();

% 
% plot(smootherTailCurveMeasure); hold on;
% plot(tmp_vel_fV_intrp);

%% ID bouts with larger gaps and remove

if ~isempty(starts)

    naned_swim_frames=[];
    
for kab=1:length(starts)
    % lost in between
    if sum(isnan(tmp_vel_fB(starts(kab):ends(kab))))>10 % >~5ms
    nanSwimFrames=ends(kab)-starts(kab);
    starts(kab)=NaN;
    ends(kab)=NaN;
    
    naned_swim_frames=[naned_swim_frames nanSwimFrames];   
    end
    
    if ~isnan(starts(kab))
    % lost just before
    if sum(isnan(tmp_vel_fB((starts(kab)-20):starts(kab))))>10 % >~28ms
    nanSwimFrames=ends(kab)-starts(kab);
    starts(kab)=NaN;
    ends(kab)=NaN;   
    naned_swim_frames=[naned_swim_frames nanSwimFrames];
    
    
    %lost just after
    elseif sum(isnan(tmp_vel_fB((ends(kab):(starts(kab)+20)))))>10 % >~28ms
    nanSwimFrames=ends(kab)-starts(kab);
    starts(kab)=NaN;
    ends(kab)=NaN;
    naned_swim_frames=[naned_swim_frames nanSwimFrames];

    
    % remove small noisy blurps of activity <50ms
    elseif length(tmp_vel_fB(starts(kab):ends(kab)))<35 % <~50ms
    nanSwimFrames=ends(kab)-starts(kab);
    starts(kab)=NaN;
    ends(kab)=NaN;
    naned_swim_frames=[naned_swim_frames nanSwimFrames];
    
    
        
    % remove other noisy blurps by using velocity thresholding
    elseif nanmax(tmp_vel_fV_intrp(starts(kab):ends(kab)))<=1.5 %1mm/s
    nanSwimFrames=ends(kab)-starts(kab);
    starts(kab)=NaN;
    ends(kab)=NaN;
    naned_swim_frames=[naned_swim_frames nanSwimFrames];
    
%     % remove other noisy blurps by using velocity thresholding
%     elseif nanmean(tmp_vel_fV_intrp(starts(kab):ends(kab)))<1 %1mm/s
%     nanSwimFrames=ends(kab)-starts(kab);
%     starts(kab)=NaN;
%     ends(kab)=NaN;
%     naned_swim_frames=[naned_swim_frames nanSwimFrames];
    
    
    end
    end
    
    
end

all_naned_swim_frames=sum(naned_swim_frames);
time_d=round(all_naned_swim_frames.*1.4286,2);

idx_NaNbouts=isnan(starts);
starts(idx_NaNbouts)=[];
ends(idx_NaNbouts)=[];

end

%%
if length(starts)>=1
    
%% plot and check filtered traces

plot(tmp_vel_unfilt); hold on;
plot(tmp_vel_fV); hold on;
plot(rad2deg(forTBF_intrp));
vline(starts,'g');
vline(ends,'r');
hold off;         

dcmObj = datacursormode;  % Turn on data cursors and return the
                          %   data cursor mode object
set(dcmObj, 'UpdateFcn', @updateFcn);  % Set the data cursor mode object update
                                       % function so it uses updateFcn.m
                                       
disp("verify traces and proceed");
%uiwait();

% %% plot representative img
% xxx=(((1:length(tmp_vel_fV)).*1.43)./1000);
% yyaxis left
% plot(xxx,tmp_vel_fV);
% yyaxis right
% plot(xxx,rad2deg(forTBF_intrp));

%% swim-burst based
for sab=1:length(starts)
    
    % bout start
    swimBouts(sab,1)=starts(sab);
    
    % bout end
    swimBouts(sab,2)=ends(sab);
        
    % bout dist in mm
    swimBouts(sab,3)=nansum(tmp_dist_fV(swimBouts(sab,1):swimBouts(sab,2)));
    
    % bout dur in ms
    swimBouts(sab,4)= (swimBouts(sab,2)-swimBouts(sab,1)).*1.4286;
    
    % mean velo
    swimBouts(sab,5)= nanmean(tmp_vel_fV_intrp(swimBouts(sab,1):swimBouts(sab,2)));
    
    % max velo
    swimBouts(sab,6)= nanmax(tmp_vel_fV_intrp(swimBouts(sab,1):swimBouts(sab,2)));
    
    % accn
    swimBouts(sab,7)= NaN;
    
    % max tail angle
    swimBouts(sab,8)= rad2deg(nanmax(abs((forMaxAngle_intrp(swimBouts(sab,1):swimBouts(sab,2))))));
      
    % tail beat frequency
    swimBouts(sab,9)= (nansum(islocalmin(forTBF_intrp(swimBouts(sab,1):swimBouts(sab,2)),'MinSeparation',10,'MinProminence',0.07)))/(swimBouts(sab,4)/1000);
    
    % max head yaw
    swimBouts(sab,10)= rad2deg(nanmax(abs(headyawFilt_intrp((swimBouts(sab,1):swimBouts(sab,2))))));
    
    % mean head yaw
    swimBouts(sab,11)= rad2deg(nanmean(abs(headyawFilt_intrp((swimBouts(sab,1):swimBouts(sab,2))))));
    
    % laterality -/+
    swimBouts(sab,12)= rad2deg(sum(headyawFilt_intrp((swimBouts(sab,1):swimBouts(sab,2)))));
    
    % angle in the begining of the bout
    swimBouts(sab,13)= rad2deg(headyawFilt_intrp((swimBouts(sab,1))));
    
    % angle at the end of the bout
    swimBouts(sab,14)= rad2deg(headyawFilt_intrp((swimBouts(sab,2))));
    
    % angular velocity
    swimBouts(sab,15)= rad2deg(sum(abs(headyawFilt_intrp((swimBouts(sab,1):swimBouts(sab,2))))))./swimBouts(sab,4);
    
   
end %end kinematics on all bouts

lskj=1;
for stg=2:(length(starts))
if sum(isnan(forTBF_intrp(ends(stg-1):starts(stg))))==0
interBoutMS(lskj)=sum(starts(stg)-ends(stg-1)).*1.4286;
lskj=lskj+1;
end
end


%% 191123 - added - 2/40K artifacts detected - negligible but still corrected

idx_fakeBouts=[];
%remove small events (<100ms) with zero displacement; to be consistent, this is done in both the fish
for sab=1:length(starts)
    if (swimBouts(sab,3)==0) & (swimBouts(sab,4)<=100)
        idx_fakeBouts=[idx_fakeBouts sab];
    end
end
swimBouts([idx_fakeBouts],:)=[];

%% operate on bout regions and ID beats - first verify IDs
% precise onset and offset of tail beat start and end is not very
% consistent across swim events and between the two species.
% Hence, i am going to use peak-to-peak half-beats for my analysis to be
% consistently reliable in identifying the half beats.
%

% ID all active regions
bb=[];
for ik=1:length(starts)
allBouts(ik).range=starts(ik):ends(ik);
bb=[bb allBouts(ik).range];
end
allActiveRange=bb';
allInactiveRange=setdiff(1:length(forTBF_intrp),allActiveRange);


% ID all half beats peaks
HBTCurveMeasure = conv(forTBF_intrp,boxcarf(3),'same'); %bcFilt=5;
HBTCurveMeasure(allInactiveRange)=NaN;
absHBTCurveMeasure=abs(HBTCurveMeasure);
%negAbsHBTcurveMeasure=absHBTCurveMeasure.*-1;
%plot(absHBTCurveMeasure);hold on; plot(forTBF_intrp); %verify filtered
%trace

[p,loc_min]=findpeaks((absHBTCurveMeasure),'MinPeakProminence',0.05,'MinPeakDistance',2);

%tmp_HBT_idx=setdiff(absHBTCurveMeasure,loc_min);
%allBeats= sort([loc_min;starts';ends']); % to include start and end from
%bout identification

allBeats=loc_min;
plt_HBT=zeros(1,length(HBTCurveMeasure))';
plt_HBT(allBeats)=1;
plot(absHBTCurveMeasure); hold on; plot(plt_HBT); hold off; %uiwait();


%% this is where half-beat based kinematics are calculated on each of the bouts - now operate on verified IDs

    HBT_dur=[];
    HBT_tbf=[];
    HBT_dist=[];
    HBT_meanVelo=[];
    HBT_maxVelo=[];
    HBT_accn=[]; 
    HBT_maxTailAngle=[]; 
    HBT_yawMax=[];
    HBT_lat=[];
    HBT_angVelo=[];
    HBT_seq=[];
    
    strtEventAll=[];
    endEventAll=[];

for jjh=1:length(starts)
    absHBTCurveMeasure_tmp=absHBTCurveMeasure(starts(jjh):ends(jjh));
    tmp_idx=starts(jjh):ends(jjh);
    
    forMaxAngle_intrp_tmp=rad2deg(forMaxAngle_intrp(starts(jjh):ends(jjh)));
    forVeloFilt_intrp_tmp=tmp_vel_fV_intrp(starts(jjh):ends(jjh));
    tmp_dist_fV_intrp_tmp=tmp_dist_fV_intrp(starts(jjh):ends(jjh));   
    forAccnFilt_tmp=forAccnFilt(starts(jjh):ends(jjh));
    headyawFilt_intrp_tmp=headyawFilt_intrp(starts(jjh):ends(jjh));

    [pks_tailPeaks,locs_tailPeaks]=findpeaks(absHBTCurveMeasure_tmp,'MinPeakProminence',0.07,'MinPeakDistance',2);
    tmp_peak_locs=tmp_idx(locs_tailPeaks);
        
    strtEvent=locs_tailPeaks(1:end-1);
    tmp_strtEvent_locs=tmp_peak_locs(1:end-1);  
    endEvent=locs_tailPeaks(2:end);
    tmp_endEvent_locs=tmp_peak_locs(2:end);
    
    %clean some events based on velo
    strtEvent(nanmax(forVeloFilt_intrp_tmp(strtEvent:endEvent))<=1.5)=[];
    tmp_strtEvent_locs(nanmax(forVeloFilt_intrp_tmp(strtEvent:endEvent))<=1.5)=[];
    endEvent(nanmax(forVeloFilt_intrp_tmp(strtEvent:endEvent))<=1.5)=[];
    tmp_endEvent_locs(nanmax(forVeloFilt_intrp_tmp(strtEvent:endEvent))<=1.5)=[];

    
    
    if length(strtEvent)>=2
        
    if rem(jjh,50)==0
    plot(abs(forTBF_intrp(starts(jjh):ends(jjh)))); hold on; plot(tmp_vel_fV_intrp(starts(jjh):ends(jjh)));vline(strtEvent,'b'),vline(endEvent,'b'); hold off;    
    disp("verify the swim burst for peaks");
    %uiwait();
    end

    for sdf=1:length(strtEvent)
        
    % dur in ms
    hbt_dur(sdf)=(endEvent(sdf)-strtEvent(sdf)).*1.4286;
    
    % HTBF
    hbt_tbf(sdf)=1000./hbt_dur(sdf);
    
    % dist in mm
    hbt_dist(sdf)=nansum(tmp_dist_fV_intrp_tmp(strtEvent(sdf):(endEvent(sdf))));

    % mean velo
    hbt_meanVelo(sdf)=nanmean(forVeloFilt_intrp_tmp(strtEvent(sdf):(endEvent(sdf))));

    % max velo
    hbt_maxVelo(sdf)=nanmax(forVeloFilt_intrp_tmp(strtEvent(sdf):(endEvent(sdf))));
    
    % mean accn
    hbt_accn(sdf)=nanmean(forAccnFilt_tmp(strtEvent(sdf):(endEvent(sdf))));
    
    % max tail angle
    hbt_maxTailAngle(sdf)=nanmax(forMaxAngle_intrp_tmp(strtEvent(sdf):(endEvent(sdf)))); % 6 frames less to only catch the peak-to-peak change and avoid the peak of the next swim cylce
 
    % max head yaw
    hbt_yawMax(sdf)=rad2deg(nanmax(headyawFilt_intrp_tmp(strtEvent(sdf):(endEvent(sdf)))));
      
    % mean head yaw - laterality
    hbt_lat(sdf)=rad2deg(nanmean(abs(headyawFilt_intrp_tmp(strtEvent(sdf):(endEvent(sdf))))));
    
    % angular velo
    hbt_angVelo(sdf)=(rad2deg(sum(abs(headyawFilt_intrp_tmp(strtEvent(sdf):(endEvent(sdf))))))./hbt_dur(sdf)).*1000; %per sec
 
    % beat # (seq of beat)
    hbt_seq(sdf)=(sdf./length(strtEvent)).*100;

    strtEventAll=[strtEventAll tmp_strtEvent_locs];
    endEventAll=[endEventAll tmp_endEvent_locs];
    
    
    if (sdf==1 & hbt_dur(sdf)>=54) | (hbt_maxVelo(sdf)<=1.5) % >36frames for the first half beat; 
    %to take care of unknown errors in detecting the first beat properly
        
    hbt_dur(sdf)=NaN;
    hbt_tbf(sdf)=NaN;
    hbt_dist(sdf)=NaN;
    hbt_meanVelo(sdf)=NaN;
    hbt_maxVelo(sdf)=NaN;
    hbt_accn(sdf)=NaN;
    hbt_maxTailAngle(sdf)=NaN;
    hbt_yawMax(sdf)=NaN;
    hbt_lat(sdf)=NaN;
    hbt_angVelo(sdf)=NaN;
    hbt_seq(sdf)=NaN;
        
    end
    
    if hbt_maxVelo(sdf)<=2 & hbt_maxVelo(sdf)>1.5
    plot(abs(forTBF_intrp(starts(jjh):ends(jjh)))); hold on; plot(forVeloFilt_intrp_tmp); vline(strtEvent(sdf),'g');vline(endEvent(sdf),'r');hold off;
    %uiwait();
    end
    
    end %end of all beats per swim event
    
    hbt_dur=hbt_dur';
    hbt_tbf=hbt_tbf';
    hbt_dist=hbt_dist';
    hbt_meanVelo=hbt_meanVelo';
    hbt_maxVelo=hbt_maxVelo';
    hbt_accn=hbt_accn';
    hbt_maxTailAngle=hbt_maxTailAngle';
    hbt_yawMax=hbt_yawMax';
    hbt_lat=hbt_lat';
    hbt_angVelo=hbt_angVelo';
    hbt_seq=hbt_seq';
    
    HBT_dur=[HBT_dur; hbt_dur];
    HBT_tbf=[HBT_tbf; hbt_tbf];
    HBT_dist=[HBT_dist; hbt_dist];
    HBT_meanVelo=[HBT_meanVelo; hbt_meanVelo];
    HBT_maxVelo=[HBT_maxVelo; hbt_maxVelo];
    HBT_accn=[HBT_accn; hbt_accn]; 
    HBT_maxTailAngle=[HBT_maxTailAngle; hbt_maxTailAngle]; 
    HBT_yawMax=[HBT_yawMax; hbt_yawMax];
    HBT_lat=[HBT_lat; hbt_lat];
    HBT_angVelo=[HBT_angVelo; hbt_angVelo];
    HBT_seq=[HBT_seq; hbt_seq];
    
    
    clearvars forTBF_intrp_tmp forMaxAngle_intrp_tmp forVeloFilt_intrp_tmp...
        tmp_dist_fV_intrp_tmp forAccnFilt_tmp headyawFilt_intrp_tmp lengthTB_ms lengthTB_frames strtEvent...
        endEvent hbt_dur hbt_tbf hbt_dist hbt_meanVelo hbt_maxVelo hbt_accn hbt_maxTailAngle hbt_yawMax hbt_lat hbt_angVelo hbt_seq;    

    end % analyse bout if atleast 2 half beats
    
end %end hbt kinematics on all bouts


%% check the peaks
%     plot(forTBF_intrp); hold on;
%     vline(strtEventAll,'g');
%     vline(endEventAll,'g'); hold off;
%     %uiwait();


%%
actTimeMS=((length(tmp_vel_fB_intrp))-(length(find(isnan(tmp_vel_fB_intrp)))+all_naned_swim_frames))*1.4286;
percentSwimmingTime=(sum(HBT_dur)/actTimeMS)*100;
swimDistancePerMin=sum(HBT_dist)/(actTimeMS/60000);


%%

swimZ_boutNorm=zscore(swimBouts);

swimZ_HBT=[HBT_dur, HBT_tbf, HBT_dist, HBT_meanVelo, HBT_maxVelo, HBT_accn, HBT_maxTailAngle, HBT_yawMax, HBT_lat, HBT_angVelo, HBT_seq];
swimZ_HBTNorm=zscore(swimZ_HBT);

%%
freeSwim(ff).species = fish_species;
freeSwim(ff).age = dpf;
freeSwim(ff).animalNum = animal_number;
freeSwim(ff).camScale = camscale_px_per_mm;
freeSwim(ff).trialLength = duration;
freeSwim(ff).rawData = tmp_data;
freeSwim(ff).totalHBTs = size(swimZ_HBT,1);
freeSwim(ff).hbtData = swimZ_HBT;
freeSwim(ff).hbtDataNorm=swimZ_HBTNorm;
freeSwim(ff).boutData=swimBouts;
freeSwim(ff).totalBouts = size(swimBouts,1);
freeSwim(ff).boutDataNorm=swimZ_boutNorm;
freeSwim(ff).percentSwimmingTime = percentSwimmingTime;
freeSwim(ff).swimDistancePerMin = swimDistancePerMin;
freeSwim(ff).actTimeMS=actTimeMS;
freeSwim(ff).interBoutMS=interBoutMS;

end

clearvars -except myFiles filePattern camscale_px_per_mm datarate_Hz freeSwim ff myDir TL

end % end file

toc
