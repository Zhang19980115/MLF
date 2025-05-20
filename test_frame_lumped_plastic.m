
%% Initialization
ModelBasic.dimension=2;
ModelBasic.nodaldof=3;
g=9.81;
Node = [];
Constrain = struct('Fix',[],'Equal',[]);
Ele = struct('Truss',[],'EBbeamcolumn',[],'Inerter',[],'Gap',[],'Cable',[],'NodeLink',[],'TSbeamcolumn',[]);
Damp = struct('rayleigh',[],'modal',[],'dashpot',[]);
Load = struct('nodstat',[],'noddyn',[],'GM',[],'TS',[],'vel',[]);


%% define variables
%unit: kN,m,ton,s

NumFloor = 6; 
H = 3.96;
B = 9.15;
Floormass = [957 957 957 957 957 1040]; 
%% create nodes 
%xcoord ycoord,zcoord, with order refecting Numbering   
for i = 1:NumFloor+1
    Node((5*i-4):(5*i),1) = [0;B;2*B;3*B;4*B];
    Node((5*i-4):(5*i),2) = (i-1)*H*ones(1,5);
end 

%% define nodal mass 
%[x y zz], with order refecting Node No.
Mass = zeros(length(Node),3);
for i = 1:NumFloor
     idx1 = (5*i+1):(5*(i+1));  %curtrent storey nodes 
     nodemass = Floormass(i)/4;
     for j = 2:4
     Mass(idx1(j),:) = [nodemass 0 0];
     end 
     for j = [1 5]
     Mass(idx1(j),:) = [nodemass/2 0 0];
     end 
end 


%% create eles 
EA_col = [0.0488 0.059 0.059 0.0488 0.012903]*2e8; %E*A
EI_col = [0.1416 0.1801 0.1801 0.1416 0.00005]*2e8; %E*I_xx
fx_col = [0.0488 0.059 0.059 0.0488 0.012903]*345000; %A*fy
fm_col = [0.07985 0.0988 0.0988 0.07985 0.000605]*345000; %z*fp    z = plastic modulus

EA_beam = [0.022387 0.022064 0.012968 0.012968 0.012968 0.012968]*2e8;
EI_beam = [0.002456 0.002052 0.000762 0.000762 0.000762 0.000762]*2e8;
fx_beam = [0.022387 0.022064 0.012968 0.012968 0.012968 0.012968]*248000;
fm_beam = [0.006801 0.006194 0.002901 0.002901 0.002901 0.002901]*158000;

EA_truss = 0.008387*2e8;
fx_truss = 0.008387*248000;

%EBbeam = #eletag #inode #jnode #Geotag  #EA #EI #fx #fm 
tag = 1;
for i = 1:NumFloor
    for j = 1:5
        idx1 = 5*(i-1)+j;
        idx2 = 5*i+j;
         Ele.EBbeamcolumn(tag,:) = [tag,idx1,idx2,1,EA_col(j),EI_col(j),fx_col(j),fm_col(j)]; %columns
        tag=tag+1;
    end 
    idx1 = (5*i+1):(5*(i+1));  %curtrent storey nodes 
    idx2 = (5*i-4):(5*i);   %lower storey nodes
    Ele.EBbeamcolumn(tag,:) = [tag,idx1(1),idx1(2),1,EA_beam(i),EI_beam(i),fx_beam(i),fm_beam(i)]; %beams
    Ele.EBbeamcolumn(tag+1,:) = [tag+1,idx1(2),idx1(3),1,EA_beam(i),EI_beam(i),fx_beam(i),fm_beam(i)]; %beams
    Ele.EBbeamcolumn(tag+2,:) = [tag+2,idx1(3),idx1(4),1,EA_beam(i),EI_beam(i),fx_beam(i),fm_beam(i)]; %beams
    Ele.EBbeamcolumn(tag+3,:) = [tag+3,idx1(4),idx1(5),1,EA_beam(i),EI_beam(i),fx_beam(i),fm_beam(i)]; %beams

    tag = tag + 4; 
end 

%Inerter = #eletag #inode #jnode #Geotag #ki #fy #fr #direction #clutch #b #cr
Ele.Inerter(1,:) = [tag,3,8,0,inf,inf,0.01,1,1,2000,500];
Ele.Inerter(2,:) = [tag+1,3,8,0,inf,inf,0.01,1,-1,2000,500];
Ele.Inerter(3,:) = [tag+2,3,8,0,inf,inf,0.01,1,0,2000,500];


%% define constrains 
%define 1 or 0 for fix or free, No. of rows decide the tag of Nodes 
Constrain.Fix = zeros(length(Node),3);
Constrain.Fix(1:5,:) = 1; 
%Constrain.Fix(6,3) = 1; 
% #master #slave #dof
% constrain.Equal(1,:)=[2 3, 1 0 0];

%% define loads 
%static nodal load 
%#node #[1xn numeric] n=NumDOF 
 % load.nodstat(5,:)=[0 0 0];  

%dynamic nodal load 
%[mx2] matrix, first column for time, second column for load magnitude
Record = '00136L';
Record_source = 'C:\Users\Yixuan\OneDrive - Imperial College London\Desktop\MLF_code_v1.1\GM_PEER\'; %%% path!!!!!!!!
Load.TS{1} = GMRecord_Shaper(Record,Record_source); 

%#node #direction=1,2,3 #TimeSeriesTag  #magnitude
% load.noddyn(1,:)=[1,1,1,0.1*g];  
Load.noddyn = [];

%#direction   #TimeSeriesTag  #magnification
Load.GM = [-1,1,0.6*g]; 


%% define damp 
%[nodes] C = alphaM*M+alphaK*K
%damp.rayleigh(1,:)=[alphaM1 alphaK1 group1];

%#nodei  #nodej  #Geotag #direction #damping coefficient
%damp.dashpot(1,:) = [5,6,1,0,1];
Damp.modal(1,:)=0.01;
Damp.modal(2,:)=0;
% damp.modal(3,:)=0.005;
% damp.modal(4,:)=0;



%% perform analysis
Optdynamic.GeoN = 1;
Optdynamic.MatN = 1;
Optdynamic.dt = 0.001;
Optdynamic.display = 1; 
Optdynamic.numstep = floor(15/Optdynamic.dt);

[out] = MLF_dynamic_solver(ModelBasic,Node,Ele,Constrain,Load,Mass,Damp,Optdynamic);


%% figure 
%calculate drift
idx = 15*(1:NumFloor)+1;
drift = 0*out.dis(:,idx);
for i = 1:length(idx)
    if i ==1
    drift(:,i) = out.dis(:,idx(1));
    else
        drift(:,i) = out.dis(:,idx(i))- out.dis(:,idx(i-1));
    end 
end 
% 
%  %close all
% 
% time-history
figure 
subplot(4,1,1)
plot(out.time,out.load.GMinput)
title('ground motion input')
grid on  
subplot(4,1,2)
plot(out.time,out.acc(:,idx)-out.load.GMinput);
grid on  
title('absolute acc')
legend ('1st storey','2nd storey')
subplot(4,1,3)
plot(out.time,out.dis(:,idx));
grid on  
title('relative disp')
subplot(4,1,4)
plot(out.time,drift);
grid on  
title('inter-storey drift')
legend('1','2','3','4','5')
% 

%% elasto-plastic behaviour 
beamtag = 6; 
forcetag = out.elelist.forcetag{beamtag,:};
deform = (out.elelist.A{beamtag,:}*out.f(:,forcetag)')'+out.defp(:,forcetag);
force = out.f(:,forcetag);

% figure
% for i = 1:length(forcetag)
%     subplot(2,2,i)
%     plot(deform(:,i),force(:,i))
%     xlabel(['strain' num2str(i)])
%     ylabel(['stress' num2str(i)])
% end 

%hysteresis history Trajectory
figureSize = [500 500 1600 500];
figure('Position', figureSize);
limx = 1.5*max(abs(deform(:,1:length(forcetag))))+1e-16;
limy = 1.5*max(abs(force(:,1:length(forcetag))))+1e-16;
speed = 10;
for i = 2:speed:(Optdynamic.numstep-speed)
    for j = 1:length(forcetag)
    subplot(1,3,j)

    if j == 3
        plot(-deform(1:speed:(i+speed),j),-force(1:speed:(i+speed),j),'Color', [0 0 0],'LineWidth',1)
        hold on 
        scatter(-deform((i+speed),j),-force((i+speed),j),'red','filled')
    else 
        plot(deform(1:speed:(i+speed),j),force(1:speed:(i+speed),j),'Color', [0 0 0],'LineWidth',1)
        hold on 
        scatter(deform((i+speed),j),force((i+speed),j),'red','filled')
    end 
    grid on 
    xlim([-limx(j) limx(j)])
    ylim([-limy(j) limy(j)])
    xlabel(['d' num2str(j)])
    ylabel(['f' num2str(j)])
     hold off
    end
        pause(0.0001)
end 


%% Clutched Inerter behaviour 
forcetag = [163 164];
force = out.f(:,forcetag);
ds = out.dis(:,22);
vs = out.vel(:,22);
vr = out.vr(:,forcetag);

%hysteresis history Trajectory
figure 
% limx = 1.5*max(abs(deform(:,1:length(forcetag))))+1e-16;
% limy = 1.5*max(abs(vs(:,1:length(forcetag))))+1e-16;
limx = out.time(end);
limy = max(abs(vs));
speed = 15;
for i = 2:speed:(Optdynamic.numstep-speed)
    plot(out.time(1:speed:(i+speed)),vs(1:speed:(i+speed)),'Color', [0.2 0.2 0.2],'LineWidth',1)
    hold on 
    plot(out.time(1:speed:(i+speed)),vr(1:speed:(i+speed),1),'r','LineWidth',1)
    plot(out.time(1:speed:(i+speed)),vr(1:speed:(i+speed),2),'b','LineWidth',1)

    scatter(out.time((i+speed)),vs((i+speed)),'k','filled')
    scatter(out.time((i+speed)),vr((i+speed),1),'r','filled')
    scatter(out.time((i+speed)),vr((i+speed),2),'b','filled')
    grid on 
    legend('Terminal Relative Speed','Flywheel A','Flywheel B')
    xlim([0 limx])
    ylim([-limy limy])
    xlabel('time(s)')
    ylabel('velocity(m/s)')
     hold off

        pause(0.0001)
end 

figure 
limx = out.time(end);
limy = max(max(abs(force)));
speed = 15;
for i = 2:speed:(Optdynamic.numstep-speed)
    plot(out.time(1:speed:(i+speed)),force(1:speed:(i+speed),1),'r','LineWidth',1)
    hold on
    plot(out.time(1:speed:(i+speed)),force(1:speed:(i+speed),2),'b','LineWidth',1)
    % scatter(out.time((i+speed)),force((i+speed),1),'r','filled')
    % scatter(out.time((i+speed)),force((i+speed),2),'b','filled')
    grid on
    xlim([0 limx])
    ylim([-limy limy])
    xlabel('time(s)')
    ylabel('Clutch force')
    hold off
     pause(0.0001)
end 


figure 
limx = max(abs(ds));
limy = max(max(abs(force)));
speed = 15;
for i = 2:speed:(Optdynamic.numstep-speed)
    plot(ds(1:speed:(i+speed)),force(1:speed:(i+speed),1),'r','LineWidth',1)
    hold on
    plot(ds(1:speed:(i+speed)),force(1:speed:(i+speed),2),'b','LineWidth',1)
    scatter(ds((i+speed)),force((i+speed),1),'r','filled')
    scatter(ds((i+speed)),force((i+speed),2),'b','filled')
    grid on
    xlim([-limx limx])
    ylim([-limy limy])
    xlabel('time(s)')
    ylabel('Clutch force')
    hold off
     pause(0.0001)
end 

figureSize =[500 500 1000 400];
figure('Position',figureSize)
limx = out.time(end);
limy = max(abs(out.load.GMinput))/9.81;
speed = 5;
for i = 2:speed:(Optdynamic.numstep-speed)
    plot(out.time,out.load.GMinput/9.81,'k','LineWidth',2)
    hold on
    plot(out.time(1:speed:(i+speed)),out.load.GMinput(1:speed:(i+speed))/9.81,'r','LineWidth',2)
    scatter(out.time((i+speed)),out.load.GMinput((i+speed))/9.81,'r','filled')
    grid on
    xlim([0 limx])
    ylim([-limy limy])
    xlabel('time(s)')
    ylabel('Ground Acceleration (g)')
    hold off
     pause(0.0001)
end 


figureSize =[500 500 1000 400];
figure('Position',figureSize)
limx = out.time(end);
limy = max(abs(out.dis(:,100)));
speed = 5;
for i = 2:speed:(Optdynamic.numstep-speed)
    plot(out.time(1:speed:(i+speed)),out.dis(1:speed:(i+speed),100)/9.81,'r','LineWidth',2)
    hold on
    scatter(out.time((i+speed)),out.dis((i+speed),100)/9.81,'r','filled')
    grid on
    xlim([0 limx])
    ylim([-limy limy])
    xlabel('time(s)')
    ylabel('Displacement')
    hold off
     pause(0.0001)
end 



 %% energy check 
idx = out.dof.dyn; 
Ek = 0.5*out.vel(:,idx).^2*diag(out.m.reduced); 
Ep = zeros(length(out.time),1); 
numele = height(out.elelist);
for i = 1:numele
A = out.elelist.A{i};
f = out.f(:,out.elelist.forcetag{i});
Ep = Ep+0.5*sum(f*A.*f,2);
end 

Es = Ek+Ep; 

figureSize = [1100 600 800 200];
figure('Position', figureSize);
plot(out.time,Es,'--')
hold on 
plot(out.time,Ek)
plot(out.time,Ep)
grid on 
legend('total energy','potential energy','kinematic energy')
title('Energy')

%%  Animation 
FigureSize = [100 100 2000 1000];
Axisrange = [-2*B 10*B;-1*H (NumFloor+6)*H];
Speed = 2.5 ;
Amplification = 10; 
Groundmove = 1;
[ploted_steps] = Animation2D(FigureSize,Axisrange,Speed,Amplification,Groundmove,out);

