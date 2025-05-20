%% input=([element_Input], [Node])
% This function will generate a [table] data for local properties of Truss element
%% output=table(#nodetag, #geotag, <chord_vector>, <alpha>, [A_local], [T_local], <mat>)

function Truss_local=ElementTruss(Truss,Node,ndim)
    if isempty(Truss) == 1
        Truss_local = table();
    else 
        nonZeroRows = any(Truss, 2);
       Truss = Truss(nonZeroRows, :);
    Numele = height(Truss);
    type = repmat({'Truss'},Numele,1);
    eletag = Truss(:,1);
    nodei = Truss(:,2);
    nodej = Truss(:,3);
    node = [nodei,nodej];
    dofi = (nodei-1)*3+(1:3);
    dofj = (nodej-1)*3+(1:3);
    geotag = Truss(:,4);

        if any([geotag == 0;geotag == 1])
        else 
        problem = num2str(find(geotag ~= 0 & geotag ~= 1));
          disp(['warning: geotag can only be 0 or 1 , cannot create tag=', problem, ' Truss element !']);
            keyboard 
        end

    coordi = Node(nodei,:);
    coordj = Node(nodej,:);    
    chord = coordj-coordi;
    Length = vecnorm(chord, 2, 2);
    alpha = atan2(chord(:,2),chord(:,1)); %calculated in assembly

    if  any(Length==0) %check input
        problem = num2str(find(Length==0));
           disp(['warning: Length=0, cannot create tag=', problem, ' Truss element !']);
            keyboard 
    end

    afpl = Length./Truss(:,5); %axial flexibility per length 
    A = (arrayfun(@(row) afpl(row), 1:height(afpl), 'UniformOutput', false))';
   
    Ti = (arrayfun(@(row) [-1,0,0], 1:height(Length), 'UniformOutput', false))';
    Tj = (arrayfun(@(row) [1,0,0], 1:height(Length), 'UniformOutput', false))';

   
    fyx = Truss(:,6);
    fy = cell(Numele,1);
    extra = cell(Numele,1);
    for i = 1:Numele
        fy{i,1} = [fyx(i) -fyx(i)];
    end    

 
    Truss_local = table(type,eletag,node,dofi,dofj,geotag,chord,alpha,A,Ti,Tj,fy,extra);
    end 
end 



