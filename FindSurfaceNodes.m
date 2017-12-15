function [outNodes] = FindSurfaceNodes(mesh, inSearchValue)
% find all nodes on a given surface
% used in SurfaceLoadSet
% class C3D8_Mesh

nodes = mesh.nodes; % mesh.nodes, all the nodes in the mesh
xnod = mesh.P(:,1); % mesh.P(:,1), the X coordinates    
searchIndices = ismember(roundNR(xnod(nodes), 4), roundNR(inSearchValue, 4));
outNodes = [nodes(searchIndices(:,1),1), nodes(searchIndices(:,2),2), ...
    nodes(searchIndices(:,3),3), nodes(searchIndices(:,4),4)...
    nodes(searchIndices(:,5),5), nodes(searchIndices(:,6),6), ...
    nodes(searchIndices(:,7),7), nodes(searchIndices(:,8),8)];

end