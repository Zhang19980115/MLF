function [ploted_steps]=Animation2D(figureSize,axisrange,speed,amplification,Groundmove,out)
t = out.time;
node = out.node;
dis = out.dis;
elelist = out.elelist;
GM = out.load.GMinput;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% !!!This function is going to plot nodes and element chords at each sampling point
%%%
%%%  <figureSize> (height x length num of pixl)
%%%  <axisrange> (height x length)
%%%  speed 
%%%  amplification 
%%%  Groundmove (0 or 1)
%%%  [node]  n x 2  [x y]  
%%%  [dis]  m x 3n   [d1 d2 d3 | d3 d4 d5 |...]
%%%  [element]  elelist
%%%  [GM]  m x 2    [t ug]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
numnode = height(node);

%coordinates of the nodes 
NodeCord=reshape(node', 1, []);

%translational dofs
TranslationalDOF=1:3*height(node);
indicesToKeep = mod(TranslationalDOF, 3) ~= 0;
TranslationalDOF = TranslationalDOF(indicesToKeep);

%displaced node coordinates
defNodeCord=NodeCord+amplification*dis(:,TranslationalDOF);

%calculate ground motion
if Groundmove==1
dt = t(2)-t(1);
ug = cumtrapz(cumtrapz(GM)*dt)*dt;

else 
    ug=zeros(length(t),2);
end 


%create figure loop
figure('Position', figureSize);
for i=[1:ceil(0.01*speed/dt):length(t), length(t)]
   
    NodeCord=defNodeCord(i,:);
    NodeCord = reshape(NodeCord, 2, numnode)';
    plot(NodeCord(:,1),NodeCord(:,2),'s','MarkerSize', 4,'LineWidth',2)
    if isempty(GM) == 0
    xlim(axisrange(1,:)+amplification*ug(i))
    ylim(axisrange(2,:));
    hold on 
    else 
    xlim(axisrange(1,:))
    ylim(axisrange(2,:));
     hold on
    end 
    
    for j=1:height(elelist)
        node1 = elelist.node(j,1);
        node2 = elelist.node(j,2);

        if strcmp(elelist.type{j},'Gap')||strcmp(elelist.type{j},'NodeLink')||strcmp(elelist.type{j},'Inerter')
            % plot(NodeCord([node1 node2],1),NodeCord([node1 node2],2),'r','LineWidth',1)%create elements chord
        elseif strcmp(elelist.type{j},'Cable')
            plot(NodeCord([node1 node2],1),NodeCord([node1 node2],2),'b','LineWidth',1)%create elements chord
        else 
            plot(NodeCord([node1 node2],1),NodeCord([node1 node2],2),'k','LineWidth',2)%create elements chord
        end 
    end 

        hold off 
   
   grid on 
   title(['t=' num2str(floor((i-1)*dt)) 's']);

 if i==1
        pause(3)
 else
    pause(0.00005)
 end 

    if i==length(t)
    else 
    cla;
    end 
end 

ploted_steps=ceil(i/ceil(0.01*speed/dt));
pause(1)
close
end 