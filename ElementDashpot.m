    
function Dashpot_local=ElementDashpot(Dashpot,Node,ndim)
    if isempty(Dashpot) == 1
        Dashpot_local = table();
    else 

    Numele = height(Dashpot);
    type = repmat({'Inerter'},Numele,1);
    eletag = Dashpot(:,1);
    nodei = Dashpot(:,2);
    nodej = Dashpot(:,3);
    node = [nodei,nodej];
    dofi = (nodei-1)*3+(1:3);
    dofj = (nodej-1)*3+(1:3);
    geotag = Dashpot(:,4);
    end 
end 