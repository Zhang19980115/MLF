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
Node(3,:)=[0.2, 0];
Node(4,:)=[0.2,1];

% define nodal mass 
 m1 = 1;
 m2 = 4;
 Mass(2,:)=[m1 0 0];
 Mass(4,:)=[m2 0 0];

%create elements
k1 = 4*pi^2 ;
k2 = 64*pi^2; 
EA = 1e5;
EI1 = k1/12;  %k = 12EI/L^3
EI2 = k2/12;
fx=inf;
fm=inf;
Ele.EBbeamcolumn(1,:) = [1,1,2,1,EA,EI1,fx,fm];
Ele.EBbeamcolumn(2,:) = [2,3,4,1,EA,EI2,fx,fm];

%Gap = #eletag #inode #jnode #Geotag #kc #fyc #direction #ug0
kc = 1e5;
fyc = inf; 
ug0 = 0.2;
Ele.Gap(1,:) = [3,2,4,0,kc,fyc,1,ug0];

% constrain 
Constrain.Fix = zeros(length(Node),3);
Constrain.Fix(1,:)=[1 1 1]; 
Constrain.Fix(3,:)=[1 1 1]; 
Constrain.Fix(:,[2 3])=1;  %lateral only 

%define loads 
Load.vel(2,:) = [-4 0 0]; %assign initial velocity
Load.vel(4,:) = [-4 0 0];

% define damping


% perform analysis
Optdynamic.GeoN = 0;
Optdynamic.MatN = 1;
Optdynamic.dt = 0.0001;
Optdynamic.display = 1;
Optdynamic.numstep = floor(8/Optdynamic.dt);

[out] = MLF_dynamic_solver(ModelBasic,Node,Ele,Constrain,Load,Mass,Damp,Optdynamic);

%% plot
idx = out.dof.dyn;

figureSize = [500 0 1000 400];
figure('Position', figureSize);

 subplot(3,1,2)
 plot(out.time,out.vel(:,idx));
 hold on 
 title('Velocity')
grid on 

subplot(3,1,3)
 plot(out.time,out.dis(:,idx));
 hold on 
 title('Displacement')
grid on 

 subplot(3,1,1)
 plot(out.time,out.f(:,7));
 hold on 
 title('Impact Force')
 grid on 
%legend('dt=0.01','dt=0.005','dt=0.001','dt=0.0005','dt=0.0001')

%% energy check 
Ek = 0.5*out.vel(:,idx).^2*[m1;m2]; 

Ep = zeros(length(out.time),1); 
numele = height(out.elelist);
for i = 1:numele
A = out.elelist.A{i};
f = out.f(:,out.elelist.forcetag{i});
Ep = Ep+0.5*sum(f*A.*f,2);
end 

Etot = Ek+Ep; 

figure 
plot(out.time,Etot,'--')
hold on 
plot(out.time,Ek)
plot(out.time,Ep)
legend('total energy','kinetic energy','potential energy')

%% Animation 
FigureSize = [500 500 500 500];
Axisrange = [-2 2;0 2];
Speed = 1 ;
Amplification = 1; 
Groundmove = 1;
[ploted_steps] = Animation2D(FigureSize,Axisrange,Speed,Amplification,Groundmove,out);

