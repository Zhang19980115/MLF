%% input=([element_Input], [Node])
% This function will generate a [table] data for local properties of Truss element
%% output=table(#nodetag, #geotag, <chord_vector>, <alpha>, [A_local], [T_local], <mat>)

function Inerter_local=ElementInerter(Inerter,Node,ndim)
    if isempty(Inerter) == 1
        Inerter_local = table();
    else 
    Numele = height(Inerter);
    type = repmat({'Inerter'},Numele,1);
    eletag = Inerter(:,1);
    nodei = Inerter(:,2);
    nodej = Inerter(:,3);
    node = [nodei,nodej];
    dofi = (nodei-1)*3+(1:3);
    dofj = (nodej-1)*3+(1:3);
    geotag = Inerter(:,4);
 
    if any([geotag == 0;geotag == 1])
    else 
        problem = num2str(find(geotag ~= 0 & geotag ~= 1));
          disp(['warning: geotag can only be 0 or 1 , cannot create tag=', problem, ' Inerter element !']);
            keyboard 
    end

    coordi = Node(nodei,:);
    coordj = Node(nodej,:);    
    chord = coordj-coordi;   
    
    ki = Inerter(:,5);
    fyx = Inerter(:,6);
    fr = Inerter(:,7);
    direction = Inerter(:,8);
    clutch = Inerter(:,9);
    b = Inerter(:,10);
    cr = Inerter(:,11);

    if  any(nodei==nodej) %check input
        problem = num2str(find(nodei==nodej)');
        disp(['warning: nodei = nodej, cannot create tag=' problem ' Inerter element!']);
            keyboard 
    end

    afpl = 1./ki; %axial flexibility per length 
    A = (arrayfun(@(row) afpl(row), 1:height(afpl), 'UniformOutput', false))';
    
    % Ti = (arrayfun(@(row) [-1,0,0], 1:Numele, 'UniformOutput', false))';
    % Tj = (arrayfun(@(row) [1,0,0], 1:Numele, 'UniformOutput', false))'; 

    fy = cell(Numele,1);
    alpha=zeros(Numele,1);
    extra = cell(Numele,1);
    Ti = cell(Numele,1);
    Tj = cell(Numele,1);
    for i = 1:Numele
       % if direction(i) == 0
       %      alpha(i) = atan2(chord(i,2),chord(i,1));  %calculated in assembly
       % else 
       %      alpha(i) = (direction(i)-1)*pi/2;
       % end  
       if direction(i) == 3
           Ti{i} = [0 0 -1];
           Tj{i} = [0 0 1];
           alpha(i) = NaN; 
       else 
           Ti{i} = [-1 0 0];
           Tj{i} = [1 0 0];
           if direction(i) == 1
                alpha(i) = (1-sign((chord(i,1)+0.00001*sign(chord(i,2)))))*pi/2;
           elseif direction(i) == 2
                alpha(i) = (1-sign((chord(i,2)+0.00001*sign(chord(i,1)))))*pi/2;
           end 
       end 


       if clutch(i) == 0
           fy{i,1} = [fyx(i) -fyx(i)]; 
       elseif clutch(i) == 1
           fy{i,1} = [fr(i)  -fyx(i)]; 
       elseif clutch(i) == -1
           fy{i,1} = [fyx(i) -fr(i) ]; 
       else 
           disp('warning: [clutch] argument can only be 1 0 or -1!');
            keyboard 
       end  
       extra{i,1} = [b(i), cr(i),clutch(i)];
    end 

    Inerter_local = table(type,eletag,node,dofi,dofj,geotag,chord,alpha,A,Ti,Tj,fy,extra);
    end 
end 




