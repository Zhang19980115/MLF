
function [out]=MLF_dynamic_solver(ModelBasic,Node,ele,constrain,load,Mass,damp,Optdynamic)
ndof = ModelBasic.nodaldof;
ndim = ModelBasic.dimension;
Numnode = height(Node);
Numdof = ndof*Numnode; 
Numstep = Optdynamic.numstep;
h = Optdynamic.dt;

%% Element 
% 【Input: ele, Node, ndim】
% 【Output: elelocal, elelist, Numele, Numforce】
% This part produce two lists of element, elelocal is partitioned based on
% type of element, elelist is the entire list. Number of element and Number
% of unknown internal forces are also derived.

%construct element list by type
elelocal.EBbeamcolumn = ElementEBbeamcolumn(ele.EBbeamcolumn,Node,ndim);
elelocal.TSbeamcolumn = ElementTSbeamcolumn(ele.TSbeamcolumn,Node,ndim);
elelocal.Truss = ElementTruss(ele.Truss,Node,ndim);
elelocal.Inerter = ElementInerter(ele.Inerter,Node,ndim);
elelocal.Gap = ElementGap(ele.Gap,Node,ndim);
elelocal.Cable = ElementCable(ele.Cable,Node,ndim);
elelocal.NodeLink = ElementNodeLink(ele.NodeLink,Node,ndim);

%construct entire element list
combined = struct2cell(elelocal);
elelist = sortrows(vertcat(combined{:}), 'eletag'); 

%check input eletag
[~, ~, idx] = unique(elelist.eletag);
counts = accumarray(idx, 1);
if any(counts>=2)
    disp('warning: repeated element tag!');
        keyboard 
end 

%create elemental force tag
Numele = height(elelist);
elelist.forcetag = cell(Numele, 1);
Numforce = 0;
for i = 1:Numele
    inc = size(elelist.A{i},1);
    elelist.forcetag{i} = (Numforce+(1:inc))';
    Numforce = Numforce + inc;
end



%% Constrain 
% 【input: constrain, ndof】
% 【output: dof.fix, dof.eql】
% This part translate input nodal constrains to related dofs
dof.fix = find(constrain.Fix'==1)';

dof.slave = [];
dof.master = [];
for i=1:height(constrain.Equal)
    slaveNode = constrain.Equal(i,2);
    masterNode = constrain.Equal(i,1);
    dofs = find(constrain.Equal(i,3:end));
    dof.slave = [dof.slave,(slaveNode-1)*ndof+dofs];
    dof.master = [dof.master,(masterNode-1)*ndof+dofs];
end 

%% Element Assemble 
T.global_cell = cell(Numele,Numnode);
for i = 1:Numele
    %if elelist.geotag(i)==1
        alpha = elelist.alpha(i);  
        col1=elelist.node(i,1);
        col2=elelist.node(i,2);
        T_alpha = [cos(alpha) sin(alpha) 0;-sin(alpha) cos(alpha) 0;0 0 1];
        T.global_cell(i,:) = {zeros(size(elelist.Ti{i}))};
        T.global_cell(i,col1) = {elelist.Ti{i}*T_alpha};
        T.global_cell(i,col2) = {elelist.Tj{i}*T_alpha};
    %end 
end 
T.global_mat = cell2mat(T.global_cell);
%apply equal constrain 
T.constrain = T.global_mat;
idx1 = dof.master;
idx2 = dof.slave;
T.constrain(:,idx1)=T.constrain(:,idx1)+T.constrain(:,idx2);  
T.constrain(:,idx2)=0;
T.constrain(:,dof.fix)=0;

%% Material behaviour 
% 【input：elelist】
% 【output: mat】
%  This part produce parameters expaining the linear (mat.A) and
%  nonlinear (mat.epp/gap/cable/inerter) behaviour of all types of elements

%local structural flexibility matrix A_local (linear part)
mat.A = blkdiag(elelist.A{:}); 
%elastic perfect plastic constrain (all forces)
mat.epp = cell2mat(elelist.fy);

%initial gap 
index = cell2mat(elelist.forcetag(strcmp(elelist.type, 'Gap')));
gap = cell2mat(elelist.extra(strcmp(elelist.type, 'Gap')));
%gap = [element force index, initial gap]
mat.gap = [index, gap];

%cable = [index, initial slackness, prestress] 
index = cell2mat(elelist.forcetag(strcmp(elelist.type, 'Cable')));
cable = cell2mat(elelist.extra(strcmp(elelist.type, 'Cable')));
mat.cable = [index, cable];

%Inerter = [index, inertance, inerter damp]
index = cell2mat(elelist.forcetag(strcmp(elelist.type, 'Inerter')));
inertance = cell2mat(elelist.extra(strcmp(elelist.type, 'Inerter')));
mat.inerter = [index, inertance];

%NodeLink = [index, yield mode] 
index = cell2mat(elelist.forcetag(strcmp(elelist.type, 'NodeLink')));
yieldmode = cell2mat(elelist.extra(strcmp(elelist.type, 'NodeLink')));
mat.link = [index, yieldmode];

%% Mass
% 【Input: Mass, Numnode, ndof】
% 【Output: dof.mass, dof.vis, mass】
% This part translate input nodal mass into related dofs and produce index
% for dofs with mass(dof.mass) and without mass(dof.vis)

%dofs first partition (mass or massless)
mass_trans = Mass';
dof.mass = find(mass_trans)';  
dof.vis = find(mass_trans==0)';

%initiation for global mass vector (all dofs)
mass.vec = zeros(Numnode*ndof,1); 

%add nodal mass to global mass matrix
mass.vec(dof.mass) = mass_trans(dof.mass)';

%modify according to equaldof (add slave dof mass to master dof mass)
mass.vec(dof.master) = mass.vec(dof.master)+mass.vec(dof.slave);
mass.vec(dof.slave) = 0;

%remove constrained dofs 
dof.mass = find(mass.vec)';  
dof.vis = find(mass.vec==0)';
idx = ~ismember(dof.vis, [dof.fix,dof.slave]);
dof.vis = dof.vis(idx); 
idx = ~ismember(dof.mass, [dof.fix,dof.slave]);
dof.mass = dof.mass(idx); 

%create mass matrix
mass.mat=diag(mass.vec);


%% Modal analysis (to obtain modal damp matrix)
% %stiffness matrix in global coordinate (no constrain)
 index = find(any(mat.A ~= 0, 1) & any(mat.A ~= 0, 2)'); 
 K.global = T.constrain(index,:)'/mat.A(index,index)*T.constrain(index,:); 
 K11 = K.global(dof.mass,dof.mass);
 K12 = K.global(dof.mass,dof.vis);
 K22 = K.global(dof.vis,dof.vis);
 K_condense = K11 - K12/K22*K12';
 M11 = mass.mat(dof.mass,dof.mass);
 [mode,omega] = eig(K_condense,M11);
 omega = sqrt(diag(omega));
 M_mode = zeros(1,length(omega));

 %sort eigen value and its vector from low to high
 [omega, sortIdx] = sort(omega);
 mode = mode(:, sortIdx);

 for i = 1:length(omega)
 M_mode(i) = mode(:,i)'*mass.mat(dof.mass,dof.mass)*mode(:,i);
 end 

 Modal.mode = mode;
 Modal.omega = omega;
 Modal.T = 2*pi./omega;
 Modal.M = M_mode;

%% Damping
% 【input: Numnode, ndof, damp, dof】
% 【output: K, mode, damp, dof.stat, dof.vis】
% This part produces damp matrix, to obtain modal damp and rayleigh
% damp, stiffness matrix K and related modes are calculated. dof.vis is
% further partition for static dofs 
c.global = zeros(Numnode*ndof); %initiation
c.rayleigh = zeros(Numnode*ndof,Numnode*ndof,height(damp.rayleigh));
c.modal = c.global;
c.dashpot = c.global;

%apply rayleigh damp 
 for i=1:height(damp.rayleigh)
    node = damp.rayleigh(i,3:end)';
    if isnan(node)
        dofs = 1:Numdof;
    else 
    dofs = reshape(((node-1)*ndof+(1:ndof))',1,[]);
    end 
    a1 = damp.rayleigh(i,1);
    a2 = damp.rayleigh(i,2);
    c.rayleigh(dofs,dofs,i) =  a1*mass.mat(dofs,dofs)+a2*K.global(dofs,dofs);
    c.global(dofs,dofs) = c.global(dofs,dofs) + c.rayleigh(dofs,dofs,i); 
 end 

%apply modal damp
idx =0;
%damp.modal(all(damp.modal == 0, 2),:)=[];
 for i=1:height(damp.modal)
     if ismember(i,1:height(omega))==0
         fprintf('defining modal damping failed, no such mode!')
         keyboard
     end 
        coefficient = damp.modal(i);
        idx = idx + 2*coefficient*omega(i)/M_mode(i)*mode(:,i)*mode(:,i)';
 end 
 c.modal(dof.mass,dof.mass) = mass.mat(dof.mass,dof.mass)*idx*mass.mat(dof.mass,dof.mass);
 c.global=c.global+c.modal;

%add dashpot
for i = 1:height(damp.dashpot)
    nodei = damp.dashpot(i,1);
    nodej = damp.dashpot(i,2);
    dofi = (nodei-1)*3+(1:3);
    dofj = (nodej-1)*3+(1:3);
    chord = Node(nodej,:)-Node(nodei,:);
    direction = damp.dashpot(i,3);
    damping = damp.dashpot(i,4);
    if direction == 0 
        alpha = atan2(chord(2),chord(1));
    else 
        alpha = (direction-1)*pi/2;
    end 
    T_alpha = [cos(alpha) sin(alpha) 0;-sin(alpha) cos(alpha) 0;0 0 1];
    T_c = [[-1 0 0]*T_alpha,[1 0 0]*T_alpha];
    c.dashpot([dofi,dofj],[dofi,dofj])= c.dashpot([dofi,dofj],[dofi,dofj]) + T_c'*damping*T_c;
end 
 c.global=c.global+c.dashpot;
 c.global(abs(c.global)<1e-10)=0;

%static dofs without dynamic behaviour (no mass, no damp)
index = find(all(c.global == 0, 1) & all(c.global == 0, 2)'); 
dof.stat = intersect(index, dof.vis); 

%remove the static dofs from dof.vis
index = ismember(dof.vis, dof.stat);
dof.vis = dof.vis(~index);
dof.dyn = sort([dof.mass,dof.vis]);

%% Load
load.dyn = zeros(Numstep+1,Numnode*ndof);  %create a matrix list of dynamic loading on all dofs
load.stat = zeros(1,Numnode*ndof); %create a vector list of static loading on all dofs
load.vel0 = zeros(1,Numnode*ndof);

%add static load 
if isempty(load.nodstat) == 0
loadvector = load.nodstat';
load.stat(:,1:length(loadvector(:))) = loadvector(:)'; 
end 



%prestress force %%%%改进迭代
if isempty(ele.Cable) == 0
P0 = T.global_mat(mat.cable(:,1),:)'*mat.cable(:,3);
load.stat = load.stat - P0';
end 

%add dynamic load 
 t = (0:h:h*Numstep)';
for i = 1:height(load.noddyn) 
    %translate input 
    load_dof = (load.noddyn(i,1)-1)*ndof+load.noddyn(i,2);
    TS = load.TS{load.GM(i,2)};
    time = TS(:,1); 
    input = TS(:,2);
    magnitude = load.noddyn(i,4);

    %interpolation
    intp_input = interp1(time,input,t,'linear',0);    
    load.dyn(:,load_dof) = intp_input*magnitude; 
end 

%add ground motion
intp_input=[]; 
for i = 1:height(load.GM)
    %translate
    gmdof = 1:ndof:Numdof;
    gmdof = dof.dyn(ismember(dof.dyn,gmdof));
    TS = load.TS{load.GM(i,2)};
    time = TS(:,1); 
    input = TS(:,2);
    
    %interpolation
    intp_input = interp1(time,input,t,'linear',0)*load.GM(i,3);    
    load.dyn(:,gmdof) = load.dyn(:,gmdof)+intp_input*mass.vec(gmdof)';
end 


if isempty(load.vel) == 0
loadvector = load.vel';
load.vel0(:,1:length(loadvector(:))) = loadvector(:)'; 
end 

load.tot = load.dyn+load.stat;
load.GMinput = intp_input;

%% MLF Dynamic Analysis 
%non-iterative parameters 
if isempty(ele.Inerter)==0
    idx_inerter = mat.inerter(:,1);
    cr = mat.inerter(:,3);
    mr = 1./mat.inerter(:,2);
else 
    idx_inerter = [];
    cr = [];
    mr = [];
end 

M = mass.mat(dof.dyn,dof.dyn);
C = c.global(dof.dyn,dof.dyn);
P1 = load.tot(:,dof.dyn)';
P2 = load.tot(:,dof.stat)';
P1 = ([zeros(length(dof.dyn),1),P1]+[P1,zeros(length(dof.dyn),1)])/2;
P1(:,[1 end])=[]; 



%non-iterative constrains 
if Optdynamic.MatN==0
    ub = inf*ones(size(mat.epp(:,1)));
    lb = -inf*ones(size(mat.epp(:,2)));
else 
    ub = mat.epp(:,1);
    lb = mat.epp(:,2);
end 

ub(idx_inerter)=mat.epp(idx_inerter,1);
lb(idx_inerter)=mat.epp(idx_inerter,2);

gamma = ones(length(mat.epp),1);
for i = idx_inerter'
    if isequal(mat.epp(i,:), [inf -inf])
        gamma(i) = 0.5;
    end 
end 
gamma = gamma(idx_inerter); 

A_leq = ones(Numforce,1);
% b_leq = inf*ones(Numforce,1);
b_leq = 2*ub;

if isempty(ele.Gap)==0
b_leq(mat.gap(:,1))=0;
end 

if isempty(ele.Cable)==0
A_leq(mat.cable(:,1))=-1;
b_leq(mat.cable(:,1))=mat.cable(:,3); %force in cable should less than invers of prestress
end 

A_leq=diag(A_leq);

%initialize lagrange multipliers 
lm.plastic = zeros(Numforce,1);
lm.eqlin = zeros(length(dof.stat),1);
lm.ineqlin = zeros(Numforce,1);

if isempty(ele.Gap)==0
lm.ineqlin(mat.gap(:,1))=mat.gap(:,2)/h; %assign initial gap
end 

if isempty(ele.Cable)==0
lm.ineqlin(mat.cable(:,1))=mat.cable(:,2)/h; %assign initial slackness
end 

%initialize iterative parameters     
B1 = T.constrain(:,dof.dyn)';
B2 = T.constrain(:,dof.stat)';
A_eq = B2;

%coefficient parameters  
y1 = (M/h+C/2)^(-1);
y2 = -0.5*y1*B1;
y3 = y1*(M/h-C/2);
y4 = gamma.*mr*h;
y5 = exp(-cr.*mr*h);
y6 = (1-gamma).*y4.*y5./gamma;
y7 = -0.5*B1'*y1;
y8 = -0.5*B1'*(eye(size(y3))+y3);
y9 =  0.5*(1+y5);
H = mat.A/h+0.25*B1'*y1*B1;
H_h = H-2*mat.A/h;
H(idx_inerter,idx_inerter) = H(idx_inerter,idx_inerter) + 0.5*diag(y4);
H_h(idx_inerter,idx_inerter) = H_h(idx_inerter,idx_inerter) + 0.5*diag(y6);
H = (H+H')/2; %avoid crash, H must be positive definite
H_h = (H_h+H_h')/2;


if Optdynamic.display == 1
 fprintf('model successfully created, start iteration!\n');
end 
 singular = 0;

%%  iteration
%initialize output variables
out.dis = zeros(Numstep+1,Numdof);
out.vel = zeros(Numstep+1,Numdof);
out.vel(1,:) = load.vel0;
out.f = zeros(Numstep+1,Numforce);
out.defp = zeros(Numstep+1,Numforce);
out.vr = zeros(Numstep+1,Numforce);
out.vr(1,idx_inerter) = out.vel(1,:)*T.constrain(idx_inerter,:)';        %match the velocity of the structure at start if initial velocity is assigned
out.defv = zeros(Numstep+1,Numforce);
out.gap = zeros(Numstep+1,Numforce);

% initialize transient output variables 
vi = out.vel(1,:)'; 
fi = out.f(1,:)';
vr_i = out.vr(1,idx_inerter)'; 

%initialize state variables
lamda.eq = zeros(Numstep+1,length(dof.stat));
lamda.eq(1,:) = -0.5*vi(dof.stat);

lamda.leq = zeros(Numstep+1,Numforce);
lamda.leq(1,:) = lm.ineqlin';

lamda.alpha = zeros(Numstep+1,Numforce);
lamda.alpha(1,:) = sqrt(b_leq-A_leq*fi)';

lamda.ub = zeros(Numstep+1,Numforce);
lamda.beta = zeros(Numstep+1,Numforce);
lamda.beta(1,:) = sqrt(ub-fi)';


lamda.lb = zeros(Numstep+1,Numforce);
lamda.gamma = zeros(Numstep+1,Numforce);
lamda.gamma(1,:) = sqrt(fi-lb)';

tic
for i = 1:Numstep
    
    %update offset <b>
    bi = H_h*fi + A_eq'*lm.eqlin - A_leq'*lm.ineqlin + lm.plastic + y7*P1(:,i) + y8*vi(dof.dyn);
    bi(idx_inerter) = bi(idx_inerter) + y9.*vr_i;

    %update equality constrain 
    b_eq = P2(:,i);
   

     %check input
      if  any([isnan(bi);abs(bi)==inf])
            singular = 1;
            info = num2str(find(isnan(bi)'+(abs(bi)==inf)'));
            fprintf(['structure collapse!!! check [b] loacation ' info ', step = ' num2str(i) '\n' ]);
            break              
      end 

        if cond(H) >= 1e20 || isnan(cond(H))
            singular = 1;
            fprintf(['matrix singular!!!, step = ' num2str(i) '\n' ]);
            break              
        end 

if Optdynamic.display == 1
      if mod(i, floor((length(t)-1) / 5)) == 0 %show process 
            fprintf('>Progress: %d%%', round(i / (length(t)-1) * 100));
      end
end 

    %options 
    tolerance=1e-12;  %play with this tolerance to balance stability, accuracy and efficiency
    algorithm = 'interior-point-convex'; 
    options = optimoptions('quadprog','Algorithm',algorithm,'Display','off','ConstraintTolerance',tolerance,'OptimalityTolerance',tolerance,'StepTolerance',tolerance); 
    options = optimoptions('quadprog','Algorithm',algorithm,'Display','off'); 

    %time-step optomization 
    [fi,~,~,~,lm] = quadprog(H,bi,A_leq,b_leq,A_eq,b_eq,lb,ub,0.95*fi,options); 
    lm.plastic = -lm.lower+lm.upper;
   
    %update ietrative parameter 
    vr_i = y6.*out.f(i,idx_inerter)' + y4.*fi(idx_inerter) + y5.*vr_i;
    vi(dof.stat) = -2*lm.eqlin;   
    vi(dof.dyn) = y1*P1(:,i)+y2*(out.f(i,:)'+fi)+y3*vi(dof.dyn);
    ui = out.dis(i,:)'+vi*h; 

    %store data in output 
    out.vel(i+1,:) = vi';
    out.f(i+1,:) = fi';
    out.dis(i+1,:) = ui';
    out.vr(i+1,idx_inerter) = vr_i; 
    out.defv(i+1,:) = 2*lm.plastic; 
    out.defp(i+1,:) =  out.defp(i,:)'+out.defv(i+1,:)'*h; 
    out.gap(i+1,:) = (h*A_leq*lm.ineqlin)';
    
    %store state variables (lambdas and slackness) 
    if isempty(lm.eqlin)==0
    lamda.eq(i,:) = lm.eqlin';
    end 
    lamda.leq(i+1,:) = lm.ineqlin';
    lamda.alpha(i+1,:) = sqrt(b_leq-A_leq*fi)';
    lamda.lb(i+1,:) = lm.lower';
    lamda.beta(i+1,:) = sqrt(ub-fi)';
    lamda.ub(i+1,:) = lm.upper';
    lamda.gamma(i+1,:) = sqrt(fi-lb)';


    %update T.constrain if GeoNon is turned on 
    if Optdynamic.GeoN == 1
        for j = 1:Numele
            if elelist.geotag(j)==1
                dofi = elelist.dofi(j,[1 2]);
                dofj = elelist.dofj(j,[1 2]);
                chord = elelist.chord(j,:) + ui(dofj)' - ui(dofi)';   
                alpha_chord = atan2(chord(:,2),chord(:,1));     
                elelist.alpha(j) = alpha_chord;  
                col1=elelist.node(j,1);
                col2=elelist.node(j,2);
                T_alpha = [cos(alpha_chord) sin(alpha_chord) 0;-sin(alpha_chord) cos(alpha_chord) 0;0 0 1];
                T.global_cell(j,col1) = {elelist.Ti{j}*T_alpha};
                T.global_cell(j,col2) = {elelist.Tj{j}*T_alpha};
            end 
        end 

        %apply equal constrain 
        T.constrain = cell2mat(T.global_cell);
        idx1 = dof.master;
        idx2 = dof.slave;
        T.constrain(:,idx1)=T.constrain(:,idx1)+T.constrain(:,idx2);
        
        %update geometric dependent variables
        B1 = T.constrain(:,dof.dyn)';
        B2 = T.constrain(:,dof.stat)';
        A_eq = B2;
        y2 = -0.5*y1*B1;
        y7 = -0.5*B1'*y1;
        y8 = -0.5*B1'*(eye(size(y3))+y3);

        H = mat.A/h+0.25*B1'*y1*B1;
        H_h = H-2*mat.A/h;
        H(idx_inerter,idx_inerter) = H(idx_inerter,idx_inerter) +  0.5*diag(y4);
        H_h(idx_inerter,idx_inerter) = H_h(idx_inerter,idx_inerter) + 0.5*diag(y6);
        H = (H+H')/2; 
        H_h = (H_h+H_h')/2;

    end 


end 
  timeused = num2str(toc);

 %% post-pocess
if singular ==1
 fprintf(['iteration failed, time used = ' timeused ' s \n' ]);
    out.f = inf(size(out.f));
    out.vel = inf(size(out.vel));
    out.dis = inf(size(out.dis));
    out.defp = inf(size(out.defp));
    out.vr = inf(size(out.vr));
    out.acc = inf(size(out.vel));
    out.defv = inf(size(out.defv));
    out.gap = inf(size(out.gap));

else 

if Optdynamic.display == 1
   fprintf(['iteration complete, time used = ' timeused ' s\n' ]);
end 
% %assign equal dof 
% out.vel(:,dof.slave) = out.vel(:,dof.master);
% out.dis(:,dof.slave) = out.dis(:,dof.master);

    if isempty(mat.cable) == 0
    out.f(:,mat.cable(:,1))=out.f(:,mat.cable(:,1))+mat.cable(:,3);
    end 

%calculate acceleration
% inter=0.5*[1,1];
% t=conv2(t,inter.','valid');
% out.f=conv2(out.f,inter.','valid');
% out.vel=conv2(out.vel,inter.','valid');
% out.dis=conv2(out.dis,inter.','valid');
% out.defp=conv2(out.defp,inter.','valid');
% out.defv = conv2(out.defv,inter.','valid');
% out.vr=conv2(out.vr,inter.','valid');
% load.GMinput=conv2(load.GMinput,inter.','valid');
out.acc = zeros(size(out.vel));
out.acc(1,:) = (out.vel(2,:)-out.vel(1,:))/h;
out.acc(2:(end-1),:) =  (out.vel(3:end,:) - out.vel(1:(end-2),:)) / (2 * h);
out.acc(end,:) = (out.vel(end,:)-out.vel((end-1),:))/h;


%collect other outputs
out.time = t; 
out.node = Node;
out.elelist = elelist;
out.dof = dof; 
out.T.global = T.global_mat;
out.T.reduced = T.constrain;
out.m.reduced = mass.mat(dof.dyn,dof.dyn);
out.m.global = mass.mat;
out.c =c; 
out.mat = mat;
out.mode = Modal;
out.load.dyn = load.dyn;
out.load.stat = load.stat;
out.load.GMinput = load.GMinput;
out.idx_inerter = idx_inerter;
end 

end 
