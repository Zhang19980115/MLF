%% input=([element_Input], [Node])
% This function will generate a [table] data for local properties of Truss element
%% output=table(#nodetag, #geotag, <chord_vector>, <alpha>, [A_local], [T_local], <mat>)

function Cable_local=ElementCable(Cable,Node,ndim)
    if isempty(Cable) == 1
        Cable_local = table();
    else 
    Numele = height(Cable);
    type = repmat({'Cable'},Numele,1);
    eletag = Cable(:,1);
    nodei = Cable(:,2);
    nodej = Cable(:,3);
    node = [nodei,nodej];
    dofi = (nodei-1)*3+(1:3);
    dofj = (nodej-1)*3+(1:3);
    geotag = Cable(:,4);

    if any([geotag == 0;geotag == 1])
    else 
        problem = num2str(find(geotag ~= 0 & geotag ~= 1));
          disp(['warning: geotag can only be 0 or 1 , cannot create tag=', problem, ' cable element !']);
            keyboard 
    end
 
    coordi = Node(nodei,:);
    coordj = Node(nodej,:);    
    chord = coordj-coordi; 
    Length = vecnorm(chord, 2, 2); 
    alpha = atan2(chord(:,2),chord(:,1)); %calculated in assembly
    
    fyt = Cable(:,6);
    ut0 = Cable(:,7);
    pt = Cable(:,8);

    if  any(Length==0) %check input
        problem = num2str(find(Length==0));
        disp(['warning: Length=0, cannot create tag=', problem, ' cable element !']);
            keyboard 
    end

    afpl = Length/Cable(:,5); %axial flexibility per length 
    A = (arrayfun(@(row) afpl(row), 1:height(afpl), 'UniformOutput', false))';
    
    Ti = (arrayfun(@(row) [-1,0,0], 1:Numele, 'UniformOutput', false))';
    Tj = (arrayfun(@(row) [1,0,0], 1:Numele, 'UniformOutput', false))'; 

    fy = cell(Numele,1);
    extra = cell(Numele,1);
    for i = 1:Numele
       fy{i,1} = [fyt(i) -fyt(i)]; 
       extra{i,1} = [ut0(i) pt(i)];
    end 

    Cable_local = table(type,eletag,node,dofi,dofj,geotag,chord,alpha,A,Ti,Tj,fy,extra);
    end 
end 


