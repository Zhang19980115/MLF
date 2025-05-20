    clc;clear;
close all

%% Initialization
ModelBasic.dimension=2;
ModelBasic.nodaldof=3;
g=9.81;
Node = [];
Constrain = struct('Fix',[],'Equal',[]);
Ele = struct('Truss',[],'EBbeamcolumn',[],'Inerter',[],'Gap',[],'Cable',[],'NodeLink',[],'TSbeamcolumn',[]);
Damp = struct('rayleigh',[],'modal',[],'dashpot',[]);
Load = struct('nodstat',[],'noddyn',[],'GM',[],'TS',[],'vel',[]);


%parameters 
H = 1; 
m = 0.2;  %concentrated mass
k = 4*pi^2; 
mb = 1; %distributed mass
N = 10;  %split
period = 2*pi*sqrt(m/k);

% create nodes 
Node(1,:)=[0, 0];   %xcoord ycoord,zcoord, with order refecting Numbering  
for i = 1:N
Node(1+i,:)=[0, i*H/N];   %xcoord ycoord,zcoord, with order refecting Numbering  
end 

%% create eles 
id = 1; 
%EBbeam = #eletag #inode #jnode #Geotag  #EA #EI #fx #fm 
EA = 1e7; 
EI = k*H^3/3;

for i = 1:N
Ele.EBbeamcolumn(i,:) = [id,i,1+i,1,EA,EI,inf,inf];
id = id +1; 
end 

%% define constrains 
%define 1 or 0 for fix or free, No. of rows decide the tag of Nodes 
Constrain.Fix = zeros(height(Node),3);
Constrain.Fix(1,:)=[1 1 1]; 

% #master #slave #dof
% constrain.Equal(1,:)=[2 3, 1 0 0];


%% define loads 
%static nodal load 
%#node #[1xn numeric] n=NumDOF 
 % load.nodstat(5,:)=[0 0 0];  

%dynamic nodal load 
%[mx2] matrix, first column for time, second column for load magnitude
Record='00136L';
Record_source='C:\Users\Yixuan\OneDrive - Imperial College London\Desktop\MLF_code_v1.1\GM_PEER\';
Load.TS{1}=GMRecord_Shaper(Record,Record_source); 

%#node #direction=1,2,3 #TimeSeriesTag  #magnitude
%load.noddyn(1,:)=[1,1,1,0.1*g];  
Load.noddyn=[];

%#direction   #TimeSeriesTag  #magnification
Load.GM=[-1,1,3*g]; 
% Load.vel(height(Node),:) = [-3 0 0]; %assign initial velocity[Node, x y zz]


%% define nodal mass 
%[x y zz], with order refecting Node No.
Mass = zeros(height(Node),3);
 Mass(end,:)=[m m 0];

 for i = 1:N
  Mass(i,:)=[mb mb 0];
 end 


%% define damp 
% alphaM1=0.01;
% alphaK1=0.00000001;
% group1=1:height(Node);
% %[nodes] C = alphaM*M+alphaK*K
% Damp.rayleigh(1,:)=[alphaM1 alphaK1 group1];

%#nodei  #nodej  #Geotag #direction #damping coefficient
%damp.dashpot(1,:) = [5,6,1,0,1];

Damp.modal(5,:)=0;
Damp.modal(9,:)=0;
Damp.modal(12,:)=0.001;
Damp.modal(14:20,:)=0.001;
%% perform analysis

% perform analysis
Optdynamic.GeoN = 1;
Optdynamic.MatN = 1;
Optdynamic.dt = 0.01;
Optdynamic.display = 1;
Optdynamic.numstep = floor(20/Optdynamic.dt);

[out] = MLF_dynamic_solver(ModelBasic,Node,Ele,Constrain,Load,Mass,Damp,Optdynamic);


%% output
 
%%  modal analysis
    figure
     hold on 
     tag = cell(1,length(out.mode.T));
     xx = 5;
for i= 1:xx
    NODE_Location = [Node,zeros(height(Node),1)]';
    NODE_Location(out.dof.mass) = NODE_Location(out.dof.mass)'+out.mode.mode(:,i);
    NODE_Location = NODE_Location([1 2],:)';
    plot(NODE_Location(:,1),NODE_Location(:,2),'o-',LineWidth=1,MarkerSize=3)
    ylim([0 1.2*H])
    xlim([-4*H 4*H])
    grid on
    tag{i} = ['mode' num2str(i) ', T=' num2str(out.mode.T(i), '%.2f') 's' ];
end 
    title('Modal Shape')
    legend(tag{1:xx})


%% figure 
idx = out.dof.dyn(end-1); 

figureSize = [500 0 500 750];
figure('Position', figureSize);

subplot(4,1,1)
 plot(t,out.load.GMinput);
 hold on 
 title('Input')
grid on 

subplot(4,1,2)
 plot(t,out.acc(:,idx));
 hold on 
 title('acc')
grid on 
subplot(4,1,3)
 plot(t,out.vel(:,idx));
 hold on 
 title('vel')
grid on 
 subplot(4,1,4)
 plot(t,out.dis(:,idx));
 hold on 
 title('dis')
grid on 
%legend('dt=0.01','dt=0.005','dt=0.001','dt=0.0005','dt=0.0001')

%% energy check 
idx = out.dof.dyn; 

Ek = 0.5*out.vel(:,idx).^2*diag(out.m.reduced); 
Ep = zeros(length(t),1); 
numele = height(elelist);
for i = 1:numele
A = elelist.A{i};
f = out.f(:,elelist.forcetag{i});
Ep = Ep+0.5*sum(f*A.*f,2);
end 

Etot = Ek+Ep; 

figureSize = [1100 600 800 200];
figure('Position', figureSize);
plot(t,Etot,'--')
hold on 
plot(t,Ek)
plot(t,Ep)
grid on 
legend('total energy','kinetic energy','potential energy')
title('Energy')

%% inerter 
% for i = 1:2
% figure
% plot(tp{i},-vel{i})
% hold on 
% plot(tp{i},vr{i})

% figure
% plot(t,-out.dis(:,4));
% hold on 
% plot(t,out.f(:,4));
% 
% 
% figure 
% plot(-out.dis(:,4),out.f(:,4));
% end 

% figure 
% for i = 2:Optdynamic.numstep
%     plot(-out.dis(i-1:i,4),out.f(i-1:i,4),'r')
%     xlim([-20 20])
%     ylim([-20 20])
%     hold on 
%      pause(0.001)
% end 





%%  ED 
% figure 
% plot(t,-out.dis(:,13),'k')
% hold on 
% plot(t,out.f(:,13))
% legend('deformation','element force')
% grid on 
% 
% % 
% figure 
% plot(-out.dis(:,13),out.f(:,13))
% 
% figure 
% subplot(2,1,1)
% plot(t,out.defu(:,4));
% hold on 
% plot(t,-out.dis(:,4));
% legend('plastic deformation','total deformation')
% 
% subplot(2,1,2)
% plot(t,-out.dis(:,4)+out.defu(:,4));
% hold on 
% plot(t,out.f(:,4))
% legend('elastic deformation','element force')

% 
% figure 
% for i = 2:Optdynamic.numstep
%     plot(out.dis(i-1:i,5),out.f(i-1:i,11),'r')
%     xlim([-0.1 0.1])
%     ylim([-2 1])
%     hold on 
%      pause(0.00000001)
% end 
% 

%% Animation 
FigureSize = [500 500 500 500];
Axisrange = [-H H;-0.5*H 1.5*H];
Speed = 1 ;
Amplification = 1; 
Groundmove = 1;
[ploted_steps] = Animation2D(FigureSize,Axisrange,Speed,Amplification,Groundmove,out);

