%% input=([element_Input], [Node])
% This function will generate a [table] data for local properties of Euler–Bernoull beam element
%% output=table(#nodetag, #geotag, <chord_vector>, alpha, [A_local], [T_local], <mat>)

function EBbeam_local=ElementEBbeamcolumn(EBbeam,Node,ndim)
    if isempty(EBbeam) == 1
        EBbeam_local = table();
    else 

    %remove zero rows
    zero_rows = all(EBbeam == 0, 2);
    EBbeam = EBbeam(~zero_rows, :);

    Numele = height(EBbeam);
    type = repmat({'EBbeamcolumn'},Numele,1);
    eletag = EBbeam(:,1);
    nodei = EBbeam(:,2);
    nodej = EBbeam(:,3);
    node = [nodei,nodej];
    dofi = (nodei-1)*3+(1:3);
    dofj = (nodej-1)*3+(1:3);
    geotag = EBbeam(:,4);

    if any([geotag == 0;geotag == 1])
    else 
        problem = num2str(find(geotag ~= 0 & geotag ~= 1));
          disp(['warning: geotag can only be 0 or 1 , cannot create tag=', problem, ' EBbeamcolumn element !']);
            keyboard 
    end

    coordi = Node(nodei,:);
    coordj = Node(nodej,:);
    chord = coordj-coordi;
    Length = vecnorm(chord, 2, 2);
    alpha = atan2(chord(:,2),chord(:,1)); %calculated in assembly
    
    if  any(Length==0) %check input
        problem = num2str(find(Length==0));
           disp(['warning: Length=0, cannot create tag=', problem, ' EBbeamcolumn element !']);
            keyboard 
    end
   
    afpl = Length./EBbeam(:,5); %axial flexibility per length 
    ffpl = Length./EBbeam(:,6); %flexural flexibility per length 
    
    A = (arrayfun(@(row) [afpl(row),0,0;0,ffpl(row)/3,-ffpl(row)/6;0,-ffpl(row)/6,ffpl(row)/3], 1:Numele, 'UniformOutput', false))';
    Ti = (arrayfun(@(row) [-1,0,0;0,1/Length(row),1;0,1/Length(row),0], 1:Numele, 'UniformOutput', false))';
    Tj = (arrayfun(@(row) [1,0,0;0,-1/Length(row),0;0,-1/Length(row),1], 1:Numele, 'UniformOutput', false))';

    fyx = EBbeam(:,7);
    fym = EBbeam(:,8);
    fy = cell(Numele,1);
    extra = cell(Numele,1);
    for i = 1:Numele
        fy{i,1} = [[fyx(i);fym(i);fym(i)],-[fyx(i);fym(i);fym(i)]];
    end 
    
   
    EBbeam_local = table(type,eletag,node,dofi,dofj,geotag,chord,alpha,A,Ti,Tj,fy,extra);

    end 
end 
