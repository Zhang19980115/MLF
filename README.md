# MLF
dynamic solver based on mixed lagrangian formalism. 

model constructing command: 
1. Node(i,:) = [x,y]  %create i_th node with coordinate x and y
2. Mass(i,:) = [mx my mz] %assign mass to node_i with component mx, my and mz
3. Constrain.Fix(i,:) = [1 0 1] %assign fix constrain to node_i on dof x y or z (1 = fix, 0 = free)
4. Constrain.Equal(1,:) = [i j, 1 0 0] %assign equaldof constrain to node_j(slave) and node_i(master) on dof x y or z(1 = equal, 0 = free)
5. Ele.EBbeamcolumn(tag,:) = [tag,node1,node2,GeoN,EA,EI,fx,fm]; %create an Euler beam element with elastic-perfect plastic behaviour (tag = unique element tag, node1 = starting node No., node2 = end node No., GeoN = 0 or 1 geometric nonlinearity, EA = axial stiffness, EI = flexural stiffness, fx = axial yield capacity, fm = flexural yield capacity)
6. Ele.Truss(tag,:) = [tag,node1,node2,GeoN,EA,fx];  %create a truss/spring element with elastic-perfect plastic behaviour (tag = unique element tag, node1 = starting node No., node2 = end node No., GeoN = 0 or 1 geometric nonlinearity, EA = axial stiffness, fx = axial yield capacity)
7. Ele.Gap(tag,:) = [tag,node1,node2,GeoN,kc,fyc,direction,ug0]; %create a frictionless 1-d contact element (tag = unique element tag, node1 = starting node No., node2 = end node No., GeoN = 0 or 1 geometric nonlinearity, kc = contact stiffness, fx = axial yield capacity, direction = 1,2 or 3 for zero-length case, ug0 = initial gap openning)
8. Ele.Cable(tag,:) = 
9. 
