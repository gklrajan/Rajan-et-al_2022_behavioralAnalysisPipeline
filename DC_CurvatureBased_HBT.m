%% 
%short timescale kinematics of DC (might be referred with its ealier name -- DT -- in
% some sections of the code)

clear; clc;

%set(0,'DefaultFigureVisible','off');

addpath(genpath('addPathToYourFunctions')); %path to the place where useful the functions are stored
cd 'addPathToWhereYourScriptsAreLocated'; %path to where the script is stored

% restoredefaultpath % This will remove any custom paths
% rehash toolboxcache
% savepath

%%  set recording settings
datarate_Hz=700;
camscale_px_per_mm=17;
TL=3;% tail length in mm

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
 %acquis_date = tmp_str{1, 1}; acquis_time = tmp_str{1, 2}; exp_type    = tmp_str{1, 3}; fish_num    = tmp_str{1, 4}; trial_num   = tmp_str{1, 6};
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
fishConfidence = tmp_data(:,8);
histogram(fishConfidence);
disp("fishConfidence histogram"); %uiwait();

idx_TE1=find(fishConfidence<150);
%tmp_data(idx_TE1,2:end)=NaN;
% can be misleading as the erraneous ID is
%often a bright spot; go to tail confidence

% based on tail confidence; sometimes a spot is accidentally detected in the absence of a fish but
% the tail confidence will still remain around zero
tailConfidence=sum(tmp_data(:,11:21),2);
histogram(tailConfidence); 
disp("tailConfidence histogram"); %uiwait();

idx_TE2=find(tailConfidence<250);
tmp_data(idx_TE2,2:end)=NaN;

tmp_frame=tmp_data(:,1);
frameNum=tmp_frame-tmp_frame(1);
frameDiff=diff(frameNum);
plot(frameDiff);
disp("freme diff - look at lost frames"); %uiwait();


%% CHECKPOINT
% CHECK FOR LOST FRAMES

% FRAME COUNTER:
frame_diff = [0; diff(tmp_data(:, 1))];

% CHECK for missing frames
idx_frame  = frameDiff > 1;   % index of missing frames
idx_lost   = find(idx_frame == 1);  % first frame in the block of missed frames


% calculate total duration of video
duration = ((length(frameNum)*1.4286)/1000)/60;

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
plot(xPos(120000:70:795000),yPos(120000:70:795000));

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


%% unfiltered raw params

dx   = [0; diff(xPos)]; % distance between two consecutive x-coordinates
dy   = [0; diff(yPos)]; % distance between two consecutive y-coordinates

tmp_dist_unfilt           = sqrt(dx.^2 + dy.^2);
tmp_dist_unfilt           = tmp_dist_unfilt./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate_Hz;         % convert to velocity in mm/s

%% prep for filtering

idx_nan     = isnan(dx);
idx_nan = find(idx_nan==1);


%% filter positions

xPos_filt = sgolayfilt(xPos,3,31);
yPos_filt = sgolayfilt(yPos,3,31);

%%
dxF   = [0; diff(xPos_filt)]; % distance between two consecutive x-coordinates
dyF   = [0; diff(yPos_filt)]; % distance between two consecutive y-coordinates

tmp_dist_filt = sqrt(dxF.^2 + dyF.^2);     % distance moved between iterations, in pixels
tmp_dist_filt = tmp_dist_filt./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_filt  = tmp_dist_filt.*datarate_Hz;  % convert to velocity in mm/s


%% velo plot - for representative slow and fast streches of swims
%dt file n=17 - 4th June 2020

% quiver(xPos_filt(100000:100:1000000),yPos_filt(100000:100:1000000),dxF(100000:100:1000000),dyF(100000:100:1000000));
% uiwait(); 

x=xPos_filt(10000:70:28700); % position per 100 ms
y=yPos_filt(10000:70:28700);

% outlier removal added on 5th june 2020 - does not change much - removes
% some high velocity noise
% approx. 0.1 to 20 becomes 0.1 to 15

%tmp_vel_filt(isoutlier(tmp_vel_filt))=NaN; 

tmp_vel_filt_gg=movmean(tmp_vel_filt,70); % 100ms moving average of velocity
v=tmp_vel_filt(10000:70:28700); % velocity per 100 ms

norm_data = (v - min(v))./(max(v)-min(v));
scatter(x, y, 30, norm_data);

nBins = 100;

% define color map
colk = hot(nBins); 

% determine color index of each point
[~, colorEdges, colorIdx] = histcounts(norm_data, nBins); 


% create figure
figure
axh = axes; 
hold(axh, 'on')
% plot data with empty markers
plot(axh, x, y, 'k-o')
% add color
scatter(x,y,35, colk(colorIdx,:), 'filled')
% add colorbar
colormap(colk)
cb = colorbar('peer', axh); 
caxis([colorEdges(1), colorEdges(end)])


%% plot and check filtered position traces
plot(xPos);
hold on;
plot(xPos_filt); disp("xpos vs filtered x-pos"); uiwait();
hold off;

plot(yPos);
hold on;
plot(yPos_filt); disp("ypos vs filtered y-pos"); %uiwait();
hold off;

plot(tmp_vel_unfilt);
hold on;
plot(tmp_vel_filt);
disp("unfilt velo vs filt velo"); uiwait();
hold off;

%%
tailSegAngles(tailSegAngles<-10)=0;
cumsumtailangles=cumsum(tailSegAngles')';

%automated selection for properly tracked segments; the others are
%discarede in teh analysis
pxVal_tail=nanmedian(tailSeg);
numsegs= 8;%length(pxVal_tail(pxVal_tail>20)); %20 now; set based on median values of few samples

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


%%
% the SNR is much lower - so we use a different approach to ID the half beats - 6 to 9
% segs


%% filter last tail seg and interpolate all signals 

rawTail=cumsumtailangles(:,numsegs);
forTBF=sgolayfilt(cumsumtailangles(:,numsegs),3,35); %8th seg - 2.4mm

forMaxAngle=abs(sgolayfilt(cumsumtailangles(:,numsegs),3,35)); %8th seg
%forVeloFilt=sgolayfilt(tmp_vel_filt,2,11);
forVeloFilt=tmp_vel_filt;
headyawFilt=sgolayfilt(headyaw_diff,2,11);

% interp small gaps ONLY - 5frames/7.1ms
forTBF_intrp=(interp1gap(forTBF,5,'linear'))';
forMaxAngle_intrp=(interp1gap(forMaxAngle,5,'linear'))';
forVeloFilt_intrp=(interp1gap(forVeloFilt,5,'linear'))';
headyawFilt_intrp=(interp1gap(headyawFilt,5,'linear'))';
tmp_dist_fV_intrp=(interp1gap(tmp_dist_filt,5,'linear'))';

% cal accn
forAccn = [0; diff(forVeloFilt_intrp)];
forAccnUnfilt = forAccn.*datarate_Hz;
forAccnFilt =sgolayfilt(forAccnUnfilt,2,15);


%%

plot(rawTail); hold on; plot(forTBF_intrp); disp("verify raw 8th segment with filtered one");
%uiwait();

%% filter the velocity trace

% filters used in the analysis
filterB = ones(1,50)./50;

idx_nan2=isnan(forVeloFilt_intrp);
forVeloFilt_intrp_a=forVeloFilt_intrp;

forVeloFilt_intrp_a(idx_nan2) = 0;
forVeloFilt_intrp_filtfilt        = filtfilt(filterB, 1, forVeloFilt_intrp_a);
forVeloFilt_intrp_filtfilt(idx_nan2) = nan; % re-insert the nan values

%plot(forVeloFilt_intrp); hold on; plot(forVeloFilt_intrp_filtfilt); plot(abs(forTBF_intrp));

%% tail angle processing from Mike Orger

threeTailSegAngles=tmp_data(1:size(tmp_data,1),27:30);%6-to-9 segs
threeTailSegAngles(threeTailSegAngles<-10)=0;
cumsumtailangles=cumsum(threeTailSegAngles')';

numsegs2=4;%length(pxVal_tail(pxVal_tail>30))

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
realSegmentAngles = fixedSegmentAngles(:,1:numsegs2);

%Filter angles
%this parameter sets the timescale at which to smooth the tail movement
bcsize=14;
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
bcFilt=14;
tailCurveMeasure = conv(superCumSegAngle(:,end),boxcarf(bcFilt),'same');

%min filter removes baseline fluctuations (needs to be well above the bout
%length)
%max filter flattens out beat variations (timescale needs to be adjusted
%according to be well below the interbout length)
smootherTailCurveMeasure = tailCurveMeasure;
%maxFilterFast(tailCurveMeasure(:,end),20) - minFilterFast(tailCurveMeasure(:,end),400);

%fix size bug of smootherTailCurveMeasure
if length(smootherTailCurveMeasure) ~= length(cumsumtailangles)
    smootherTailCurveMeasure(end) = [];
end

% filter used on the curvature measure
filterB = ones(1,200)./200;

idx_nan3=isnan(smootherTailCurveMeasure);
smootherTailCurveMeasure_a=smootherTailCurveMeasure;

smootherTailCurveMeasure_a(idx_nan3) = 0;
smootherTailCurveMeasure_a_filtfilt        = filtfilt(filterB, 1, smootherTailCurveMeasure_a);
smootherTailCurveMeasure_a_filtfilt(idx_nan3) = nan; % re-insert the nan values

%manual threshold
allbout=smootherTailCurveMeasure_a_filtfilt>0.13;
allstart=max([0 diff(allbout')],0);
allend=max([0 -diff(allbout')],0);

%% ID regions to calculate bout associated kinematics

starts=find(allstart);
ends=find(allend);

% first calculate ALL the active regions ID'd based on tail curvature
% irrespective of the length of the movements before filterting for longer
% events below
swmFrm=ends-starts;
swmTme=round(sum(swmFrm).*1.4286,2);



if starts(1)>ends(1)
ends(1)=[]; % to remove already init swim
end

if starts(end)>ends(end)
starts(end)=[]; % to remove half swim in the end
end

%
% for jk=1:length(starts)
%     chk=nanmean(forVeloFilt_intrp_filtfilt(starts(jk):ends(jk)));
%     
%     if chk<1
%         starts(jk)=NaN;
%         ends(jk)=NaN;
%     end
% end
% 
% starts(isnan(starts))=[];
% ends(isnan(ends))=[];

% remove events smaller than 200 ms - around 4 cycles

naned_swim_frames=[];

for kjk=1:length(starts)

    if ends(kjk)-starts(kjk)<=140
    nanSwimFrames=ends(kjk)-starts(kjk);
    naned_swim_frames=[naned_swim_frames nanSwimFrames];

    starts(kjk)=NaN;
    ends(kjk)=NaN;
    end
    
end


idx_NaNbouts=isnan(starts);
starts(idx_NaNbouts)=[];
ends(idx_NaNbouts)=[];

all_naned_swim_frames=sum(naned_swim_frames);
time_d=round(all_naned_swim_frames.*1.4286,2);


plot(forTBF_intrp); hold on; plot(smootherTailCurveMeasure_a_filtfilt);
vline(starts,'g'); vline(ends,'r'); %%%%%%uiwait();


%% identify swim islands - start and end of well tracked regions
% 
% % Get logical map of where the nans are.
% trackedLocations = ~isnan(forTBF_intrp);
% labeledSwims = bwlabel(trackedLocations);
% 
% % Find their areas
% measurements = regionprops(labeledSwims, 'Area', 'PixelIdxList');
% allAreas = [measurements.Area];
% 
% % Find label of those areas 140 frames/ 200ms or more / 4*cycle
% bigAreas = find(allAreas>140);
% 
% % find starts and ends of well tracked swim blocks
% for trk = 1 : length(bigAreas)
%   tmp_block=measurements((bigAreas(trk))).PixelIdxList;
% 
%   % first element of the swim block
%   trackStart(trk)=tmp_block(1); %
%   
%   % last element of the swim block
%   trackEnd(trk)=tmp_block(end);
%   
% end



  %% now find swim events in the well tracked regions and cal hbt

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
    
    
    %% plot representative img

xxx=(((1:length(forVeloFilt_intrp_filtfilt)).*1.43)./1000);
yyaxis left
plot(xxx,forVeloFilt_intrp_filtfilt);
yyaxis right
plot(xxx,rad2deg(forTBF_intrp));
    
    
  if length(starts)>=1
  
    for swm = 1 : length(starts)    
    
    tmp_forTBF_intrp=forTBF_intrp(starts(swm):ends(swm));
    forMaxAngle_intrp_tmp=rad2deg(forMaxAngle_intrp(starts(swm):ends(swm)));
    forVeloFilt_intrp_tmp=forVeloFilt_intrp_filtfilt(starts(swm):ends(swm));
    tmp_dist_fV_intrp_tmp=tmp_dist_fV_intrp(starts(swm):ends(swm));   
    forAccnFilt_tmp=forAccnFilt(starts(swm):ends(swm));
    headyawFilt_intrp_tmp=headyawFilt_intrp(starts(swm):ends(swm));
    
    
    [pks_tmp,locs_tmp]=findpeaks(abs(tmp_forTBF_intrp),'MinPeakProminence',0.02,'MinPeakDistance',2);
    
    % check the filtered velocity at all the tail peaks and filtere with
    % a thresdhold velo - this will make sure that all the noisy peaks are
    % eliminated while the occasional tiny tail peaks embedded in a large
    % movement will still be retained
    
    id_falsePeaks=locs_tmp(forVeloFilt_intrp_filtfilt(locs_tmp)<1);
    locs_tailPeaks=setdiff(locs_tmp,id_falsePeaks);
        
%     plot(abs(forTBF_intrp(trackStart(swm):trackEnd(swm)))); hold on; plot(forVeloFilt_intrp_filtfilt(trackStart(swm):trackEnd(swm)));vline(locs_tailPeaks); hold off;    
%     disp("verify the tracked island for peaks");
%     uiwait();
    
    
    strtEvent=[locs_tailPeaks(1:end-1)];
    endEvent=[locs_tailPeaks(2:end)];

    lengthTB_frames=[diff(locs_tailPeaks)];
    lengthTB_ms=(((diff(locs_tailPeaks)).*1.4286));
    
    idx_swimSwitch=lengthTB_frames>=50; %70 ms between half beats
    strtEvent(idx_swimSwitch)=[];
    endEvent(idx_swimSwitch)=[];
    
    
    if length(strtEvent)>=2
    
    
    if rem(swm,10)==0
    plot(abs(forTBF_intrp(starts(swm):ends(swm)))); hold on; plot(forVeloFilt_intrp_filtfilt(starts(swm):ends(swm)));vline(strtEvent,'g'),vline(endEvent,'r'); hold off;    
    disp("verify the swim burst for peaks");
    %%%%%%uiwait();
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
    hbt_seq(sdf)=sdf;
    
    
    if (sdf==1 & hbt_dur(sdf)>=54) % >36frames for the first half beat; 
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
    
    
    end %check if there are atleast two beats in the ID'd swim event
    
    end %end of all ID'd swim events in a file
    
    
%%
actTimeMS=((length(forVeloFilt_intrp_filtfilt))-(length(find(isnan(forVeloFilt_intrp_filtfilt)))+all_naned_swim_frames))*1.4286;
percentSwimmingTime=(sum(HBT_dur)/actTimeMS)*100;
swimDistancePerMin=sum(HBT_dist)/(actTimeMS/60000);

%%
swimD_HBT=[HBT_dur, HBT_tbf, HBT_dist, HBT_meanVelo, HBT_maxVelo, HBT_accn, HBT_maxTailAngle, HBT_yawMax, HBT_lat, HBT_angVelo, HBT_seq];
swimD_HBTNorm=zscore(swimD_HBT);

%%
freeSwim(ff).species = fish_species;
freeSwim(ff).age = dpf;
freeSwim(ff).animalNum = animal_number;
freeSwim(ff).camScale = camscale_px_per_mm;
freeSwim(ff).trialLength = duration;
freeSwim(ff).rawData = tmp_data;

freeSwim(ff).totalHBTs = size(swimD_HBT,1);
freeSwim(ff).hbtData = swimD_HBT;
freeSwim(ff).hbtDataNorm=swimD_HBTNorm;

freeSwim(ff).percentSwimmingTime = percentSwimmingTime;
freeSwim(ff).swimDistancePerMin = swimDistancePerMin;
freeSwim(ff).actTimeMS=actTimeMS;

    
  end %check if ID'd swim event is long enough to accomodate beats - should always be true as we screen for this already

clearvars -except myFiles filePattern camscale_px_per_mm datarate_Hz freeSwim ff myDir TL

end % end of file

toc;
