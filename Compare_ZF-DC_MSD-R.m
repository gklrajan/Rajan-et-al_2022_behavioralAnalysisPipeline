% Analyse mean square displacement and resultant vector length for Danionella and Zebrafish trajectory.
% 

datapath='insertYourDataDirectory';

pix=17;     % pixel size in mm
frequency=700; % framerate of acquisition

downsample=70; % data are filtered then downsampled by this factor for the analysis.
tau=downsample/700; % time period between datapoint after downsampling (0.1s here)
center=[400,400];
Radius=350; % estimated center and radius of the ROI, to eliminate data too close to the wall
dt_v=2; % time interval over which to compute heading orientation (in number of tau)
delta=[0:1:200]; % time indices at which we evaluate MSD (i.e. 0 to 20s every 0.1s)


%% Start with Danionella
load([datapath 'DT.mat']);
Ndt=length(dt);

MSD=NaN(length(delta),Ndt); % will get the MSD per fish
MSDsem=NaN(length(delta),Ndt); % will get the sem of MSD per fish
R=NaN(length(delta),Ndt); % will get the mean resultant vector length per fish
Rdiff=NaN(length(delta),Ndt); % transverse part of the orientation
Rsem=NaN(length(delta),Ndt); % will get the sem of the resultant vector length per fish
Npoints_dt=zeros(Ndt,1); % number of datapoints within ROI

for n=1:Ndt; % loop on fish
 display(['working on fish ' num2str(n) ' out of ' num2str(Ndt)])
 x=dt(n).xPos;
 y=dt(n).yPos;
 x=sgolayfilt(x,2,2*downsample+1); % filter trajectories over twice the sampling window
 y=sgolayfilt(y,2,2*downsample+1);
 x=x([1:downsample:end]); % downsample
 y=y([1:downsample:end]);
 w=find(((x-center(1)).^2+(y-center(2)).^2)>Radius^2); % detect positions too close to the border
 IndIn=ones(length(x),1);
 IndIn(w)=0; % IndIn is one when the fish is inside the ROI, zero otherwise

% separate the data into a discrete collection of trajectories within the ROI 
[L,Ntraj]=bwlabel(IndIn);    % Ntraj is the number of trajectories
Ltraj=zeros(Ntraj,1);        % Ltraj gives the length of each trajectory
for traj=1:Ntraj;
    Ltraj(traj)=length(find(L==traj));
end

Xsp=NaN(Ntraj, max(Ltraj));         % x,y position split in trajectories
Ysp=NaN(Ntraj, max(Ltraj));
theta=NaN(Ntraj, max(Ltraj)-dt_v);  % heading orientation

Npoints_dt(n)=length(find(IndIn==1 & isnan(x)==0));
Npoints_dt(n)

for traj=1:Ntraj;
    Xsp(traj,1:Ltraj(traj))=x(find(L==traj));
    Ysp(traj,1:Ltraj(traj))=y(find(L==traj));
    dXsp=Xsp(traj,(dt_v+1):end)-Xsp(traj,1:end-dt_v);
    dYsp=Ysp(traj,(dt_v+1):end)-Ysp(traj,1:end-dt_v);
    w=find(dXsp.^2+dYsp.^2)>100;
    theta(traj,w)=atan2(dYsp(w),dXsp(w));
    quiver(Xsp(traj,w),Ysp(traj,w),dXsp(w),dYsp(w));
    plot(Xsp(traj,w),Ysp(traj,w),'*');
    %plot(theta(traj,w)-theta(traj,w(1)));
   
%     theta(traj,:)=atan2(dYsp,dXsp);
%     quiver(Xsp(traj,1:end-dt_v),Ysp(traj,1:end-dt_v),10*cos(theta(traj,:)),10*sin(theta(traj,:)));
      hold on
    
end

% plot trajectories
%plot(Xsp',Ysp','*-')

for k=1:length(delta);
        p=delta(k);
        dx=Xsp(:,(p+1):end)-Xsp(:,1:end-p);
        dy=Ysp(:,(p+1):end)-Ysp(:,1:end-p);
        MSD(k,n)= nanmean(dx(:).^2+dy(:).^2);
        Ndata=sum(isfinite(dx(:))==1);
        if Ndata >=1;
            semMSD(k,n)=nanstd(dx(:).^2+dy(:).^2)/sqrt(Ndata);
        end
        
        dtheta=theta(:,(p+1):end)-theta(:,1:end-p);
        % R(k,n)=sqrt(nanmean(cos(dtheta(:)))^2+nanmean(sin(dtheta(:)))^2);
        R(k,n)=nanmean(cos(dtheta(:)));
        Rdiff(k,n)=sqrt(nanmean(sin(dtheta(:)))^2);
        Ndata=sum(isfinite(dtheta(:))==1);
%         if Ndata >=1;
%             semR(k,n)=nanstd(cos(dtheta(:)))^2+nanmean(sin(dtheta(:)))^2)/sqrt(Ndata);
%         end
end

end

MSDdt=MSD/pix^2; % convert in mm^2
Rdt=R;
clear('dt'); % takes a lot of RAM


% plotting
% plot all MSD:
figure
hold off
for n=1:Ndt
    plot(delta*tau, MSDdt(:,n));
    hold on;
end
xlim([0,10]);
xlabel('time (second)');
ylabel('MSD (mm^2)');
xlim([0 10]);
set(gca,'FontSize',16);
title('Mean Square Displacement (all fish) - Danionella');

% plot all R
figure
hold off
for n=1:Ndt
    plot(delta*tau,Rdt(:,n));
    hold on
end
xlabel('time (second)');
ylabel('R');
xlim([0 10]);
set(gca,'FontSize',16);
title('Resultant Vector Length (all fish) - Danionella');

% plot mean MSD
figure
errorbar(delta*tau,nanmean(MSDdt,2), nanstd(MSDdt,1,2)/sqrt(Ndt), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('MSD (mm^2)');
xlim([0 10]);
set(gca,'FontSize',16);
title(['Mean Square Displacement (N=' num2str(Ndt) ') - Danionella']);

% plot mean R
figure
errorbar(delta*tau,nanmean(Rdt,2), nanstd(Rdt,1,2)/sqrt(Ndt), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('R');
xlim([0 10]);
set(gca,'FontSize',16);
title(['Resultant vector length (N=' num2str(Ndt) ') - Danionella']);

%% Now zebrafish
load([datapath 'ZF.mat']);
Nzf=length(zf);
figure

MSD=NaN(length(delta),Nzf); % will get the MSD per fish
R=NaN(length(delta),Nzf); % will get the mean R vector length per fish
Npoints_zf=zeros(Nzf,1); % number of datapoints within ROI

for n=1:Nzf; % loop on fish
 display(['working on fish ' num2str(n) ' out of ' num2str(Nzf)])
 x=zf(n).xPos;
 y=zf(n).yPos;
 x=sgolayfilt(x,2,2*downsample+1);
 y=sgolayfilt(y,2,2*downsample+1);
 x=x([1:downsample:end]);
 y=y([1:downsample:end]);
 w=find(((x-center(1)).^2+(y-center(2)).^2)>Radius^2); % detect of positions too close to the border
 IndIn=ones(length(x),1);
 IndIn(w)=0; % IndIn is one when the fish is inside the ROI, zero otherwise
 Npoints_zf(n)=length(find(IndIn==1 & isnan(x)==0));
 Npoints_zf(n)
 
% separate the data into a discrete collection of trajectories within the ROI 
[L,Ntraj]=bwlabel(IndIn); % Ntraj is the number of trajectories
Ltraj=zeros(Ntraj,1); % Ltraj gives the length of each trajectory
for traj=1:Ntraj;
    Ltraj(traj)=length(find(L==traj));
end

Xsp=NaN(Ntraj, max(Ltraj));
Ysp=NaN(Ntraj, max(Ltraj));
theta=NaN(Ntraj, max(Ltraj)-dt_v);  % for heading orientation

for traj=1:Ntraj;
    Xsp(traj,1:Ltraj(traj))=x(find(L==traj));
    Ysp(traj,1:Ltraj(traj))=y(find(L==traj));
    dXsp=Xsp(traj,(dt_v+1):end)-Xsp(traj,1:end-dt_v);
    dYsp=Ysp(traj,(dt_v+1):end)-Ysp(traj,1:end-dt_v);
    w=find(dXsp.^2+dYsp.^2)>100;
    theta(traj,w)=atan2(dYsp(w),dXsp(w));
    %theta(traj,:)=atan2(dYsp,dXsp);
    %quiver(Xsp(traj,1:end-dt_v),Ysp(traj,1:end-dt_v),10*cos(theta(traj,:)),10*sin(theta(traj,:)));
    %quiver(Xsp(traj,w),Ysp(traj,w),10*cos(theta(traj,w)),10*sin(theta(traj,w)));
    quiver(Xsp(traj,w),Ysp(traj,w),dXsp(w),dYsp(w));
    
    hold on
%     plot(Xsp(traj,1:end-dt_v),Ysp(traj,1:end-dt_v),'*-')
%     hold on;
%     quiver(Xsp(traj,1:end-dt_v),Ysp(traj,1:end-dt_v),10*cos(theta(traj,:)),10*sin(theta(traj,:)));
%     hold off
%     pause(0.5)
end

% plot trajectories
%plot(Xsp',Ysp','*-')
%pause(.1)

for k=1:length(delta);
        p=delta(k);
        
        dx=Xsp(:,(p+1):end)-Xsp(:,1:end-p);
        dy=Ysp(:,(p+1):end)-Ysp(:,1:end-p);
        MSD(k,n)= nanmean(dx(:).^2+dy(:).^2);
        Ndata=sum(isfinite(dx(:))==1);
        if Ndata >=1;
            stdMSD(k,n)=nanstd(dx(:).^2+dy(:).^2)/sqrt(Ndata);
        end
        
        dtheta=theta(:,(p+1):end)-theta(:,1:end-p);
        %R(k,n)=sqrt(nanmean(cos(dtheta(:)))^2+nanmean(sin(dtheta(:)))^2);
        R(k,n)=nanmean(cos(dtheta(:)));
        Ndata=sum(isfinite(dtheta(:))==1);
end

end

MSDzf=MSD/pix^2; % convert in mm^2
Rzf=R;
clear('zf'); % takes a lot of RAM

% plotting
% plot all MSD:
figure
hold off
for n=1:Nzf
    plot(delta*tau, MSDzf(:,n));
    hold on;
end
xlim([0,10]);
xlabel('time (second)');
ylabel('MSD (mm^2)');
xlim([0 10]);
set(gca,'FontSize',16);
title('Mean Square Displacement (all fish) - Zebrafish');

% plot all R
figure
hold off
for n=1:Nzf
    plot(delta*tau,Rzf(:,n));
    hold on
end
xlabel('time (second)');
ylabel('R');
xlim([0 10]);
set(gca,'FontSize',16);
title('Resultant Vector Length (all fish) - Zebrafish');

% plot mean MSD
figure
errorbar(delta*tau,nanmean(MSDzf,2), nanstd(MSDzf,1,2)/sqrt(Nzf), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('MSD (mm^2)');
xlim([0 10]);
set(gca,'FontSize',16);
title(['Mean Square Displacement (N=' num2str(Nzf) ') - Zebrafish']);


% plot mean R
figure
errorbar(delta*tau,nanmean(Rzf,2), nanstd(Rzf,1,2)/sqrt(Nzf), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('R');
xlim([0 10]);
set(gca,'FontSize',16);
title(['Resultant vector length (N=' num2str(Nzf) ') - Zebrafish']);
1;

%% Comparative plots
figure
s1=subplot(1,2,1);

% plot mean MSD
e1=errorbar(delta*tau,nanmean(MSDzf,2), nanstd(MSDzf,1,2)/sqrt(Nzf), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
hold on;
e2=errorbar(delta*tau,nanmean(MSDdt,2), nanstd(MSDdt,1,2)/sqrt(Ndt), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('MSD (mm^2)');
xlim([0 10]);
set(gca,'FontSize',16);
title('Mean Square Displacement');
set(e1,'DisplayName','ZF');
set(e2,'DisplayName','DT');
legend1 = legend(s1,'show');
set(legend1,'Position',[0.25 0.8 0 0]);

s2=subplot(1,2,2);

% plot mean R
e1=errorbar(delta*tau,nanmean(Rzf,2), nanstd(Rzf,1,2)/sqrt(Nzf), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
hold on;
e2=errorbar(delta*tau,nanmean(Rdt,2), nanstd(Rdt,1,2)/sqrt(Ndt), ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('R');
xlim([0 10]);
set(gca,'FontSize',16);
title('Resultant Vector Length');
set(e1,'DisplayName','ZF');
set(e2,'DisplayName','DT');
legend1 = legend(s2,'show');
set(legend1,'Position',[0.75 0.8 0 0]);

%%
MSDmdt=nanmean(MSDdt,2);% mean MSD for danionella
MSDmzf=nanmean(MSDzf,2);% mean MSD for zebrafish
semMSDdt=nanstd(MSDdt')/sqrt(Ndt);
semMSDzf=nanstd(MSDzf')/sqrt(Nzf);

Ddt=nanmean(sqrt(MSDdt),2); % mean travelled distance
Dzf=nanmean(sqrt(MSDzf),2); % mean travelled distance
semDdt=nanstd(Ddt')/sqrt(Ndt);
semDzf=nanstd(Dzf')/sqrt(Nzf);

Rmdt=nanmean(Rdt,2); % angular correlation
Rmzf=nanmean(Rzf,2);
semRdt=nanstd(Rdt')/sqrt(Ndt);
semRzf=nanstd(Rzf')/sqrt(Nzf);

Dbaldt=[0 cumsum(Rmdt(1:end-1)*Ddt(2))']; % ballistic component (Ddt(2) is ?the instanteneous velocity
Dbalzf=[0 cumsum(Rmzf(1:end-1)*Dzf(2))'];


figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% plot mean MSD
e1=errorbar(delta*tau,MSDmzf,semMSDzf, ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
hold on;
e2=errorbar(delta*tau,MSDmdt,semMSDdt, ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
p1=plot(delta*tau,Dbalzf.^2,'linewidth',2);
p2=plot(delta*tau,Dbaldt.^2,'linewidth',2);
xlabel('time (second)');
ylabel('MSD (mm^2)');
xlim([0 12]);
set(gca,'FontSize',16);
title('Mean Square Displacement');
set(e1,'DisplayName','zebrafish');
set(e2,'DisplayName','Danionella');
set(p1,'DisplayName', 'zebrafish ballistic component');
set(p2,'DisplayName', 'Danionella ballistic component');
legend1 = legend(axes1,'show');
xlim(axes1,[0 12]);
box(axes1,'on');
set(axes1,'FontSize',16);
set(legend1,...
    'Position',[0.25 0.6 0.199764150943396 0.29090909090909]);

figure
% plot mean R
e1=errorbar(delta*tau,Rmzf, semRzf, ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
hold on;
e2=errorbar(delta*tau,Rmdt, semRdt, ...
    'MarkerSize',10,'Marker','*','LineWidth',1);
xlabel('time (second)');
ylabel('R');
xlim([0 12]);
set(gca,'FontSize',16);
title('heading orientation correlation');
set(e1,'DisplayName','ZF');
set(e2,'DisplayName','DT');


