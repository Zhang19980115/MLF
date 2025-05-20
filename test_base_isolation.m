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

% create nodes 
Node(1,:)=[0, 0];   %xcoord ycoord,zcoord, with order refecting Numbering  
Node(2,:)=[0, 1];
Node(3,:)=[0, 11];
Node(4,:)=[1,1];
Node(5,:)=[-1,1];

mb = 0.1;
m2 = 1;
% define nodal mass 
 Mass(2,:)=[mb 0 0];
 Mass(3,:)=[m2 0 0];

%create elements
kb = 10 ;
k = 4*pi^2; 
EA = 1e5;
EI1 = kb/12;  %EI = kL^3/12
EI2 = k*10^3/12;
fx=inf;
fm=inf;
Ele.EBbeamcolumn(1,:) = [1,1,2,0,EA,EI1,fx,fm];
Ele.EBbeamcolumn(2,:) = [2,2,3,1,EA,EI2,fx,fm];

%Gap = #eletag #inode #jnode #Geotag #kc #fyc #direction #ug0
kc = 1e4;
fyc = inf; 
ug0 = 0.2;
Ele.Gap(1,:) = [3,4,2,0,kc,fyc,1,ug0];
Ele.Gap(2,:) = [4,2,5,0,kc,fyc,1,ug0];
% constrain 
Constrain.Fix = zeros(length(Node),3);
Constrain.Fix(1,:)=[1 1 1]; 
Constrain.Fix(4,:)=[1 1 1]; 
Constrain.Fix(5,:)=[1 1 1]; 
Constrain.Fix(:,[2 3])=1;  %lateral only 

%define loads 
Load.vel(3,:) = [-1 0 0]; %assign initial velocity

% define damping


% perform analysis
Optdynamic.GeoN = 0;
Optdynamic.MatN = 1;
Optdynamic.dt = 0.001;
Optdynamic.display = 1;
Optdynamic.numstep = floor(3/Optdynamic.dt);

[out] = MLF_dynamic_solver(ModelBasic,Node,Ele,Constrain,Load,Mass,Damp,Optdynamic);


%% plot 
figureSize = [500 100 1000 400];
figure('Position', figureSize);
idx = out.dof.dyn;
 subplot(3,1,2)
 plot(out.t,out.vel(:,idx));
 hold on 
 title('acc')
 grid on 


subplot(3,1,3)
 plot(out.t,out.dis(:,idx));
 hold on 
 title('vel')
grid on 

 subplot(3,1,1)
 plot(out.t,out.f(:,7));
 hold on 
 title('impact force')
 grid on 
%legend('dt=0.01','dt=0.005','dt=0.001','dt=0.0005','dt=0.0001')

%% energy check 
Ek = 0.5*out.vel(:,idx).^2*[mb;m2]; 
Ep = zeros(length(out.t),1); 
numele = height(out.elelist);

for i = 1:numele
A = out.elelist.A{i};
f = out.f(:,out.elelist.forcetag{i});
Ep = Ep+0.5*sum(f*A.*f,2);
end 

Etot = Ek+Ep; 

figureSize = [1200 800 500 400];
figure('Position', figureSize);
plot(out.t,Etot,'--')
hold on 
plot(out.t,Ek)
plot(out.t,Ep)
legend('total energy','kinetic energy','potential energy')
 
%% Animation 
FigureSize = [500 500 500 500];
Axisrange = [-2 2;0 15];
Speed = 1 ;
Amplification = 1; 
Groundmove = 1;
[ploted_steps] = Animation2D(FigureSize,Axisrange,Speed,Amplification,Groundmove,out);

