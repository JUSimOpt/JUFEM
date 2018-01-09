function h = vizBCandLoads(model)
scaleBC = 0.3;
scaleLoad = 1;

% TODO: input arguments
% TODO: generelize for 2D

mesh = model.mesh;

xfigure
h = mesh.vizMesh();
h.patch.FaceColor = 'w';
light

dofs = model.dofs;

hold on
P = model.mesh.P;
BC = model.mesh.BoundarieConditions;

nB = length(BC);
h = [];
for i = 1:nB
    ind = BC(i).NodeSet;
    I = ones(length(ind),1);
    if isfield(BC(i),'U1')
        h(end+1) = quiver3(P(ind,1),P(ind,2),P(ind,3),I,I*0,I*0,scaleBC,'LineWidth',3,'Color','r','ShowArrowHead','off');
    end
    if isfield(BC(i),'U2')
        h(end+1) = quiver3(P(ind,1),P(ind,2),P(ind,3),I*0,I,I*0,scaleBC,'LineWidth',3,'Color','r','ShowArrowHead','off');
    end
    if isfield(BC(i),'U3')
        h(end+1) = quiver3(P(ind,1),P(ind,2),P(ind,3),I*0,I*0,I,scaleBC,'LineWidth',3,'Color','r','ShowArrowHead','off');
    end
end


Loads = model.mesh.Loads;
nL = length(Loads);
for i = 1:nL
    nodes = Loads(i).Set;
    nele = size(nodes,1);
    if Loads(i).Type == LoadType.PointLoad
        Xm = P(nodes,:);
    else
        %        find mean point
        Xm = zeros(nele,dofs);
        for iel = 1:nele
            iv = nodes(iel,:);
            Xm(iel,:)=mean(P(iv,:),1);
        end
    end
    dir = Loads(i).Direction;
    I = ones(nele,1);
    h(end+1) = quiver3(Xm(:,1),Xm(:,2),Xm(:,3),dir(1)*I,dir(2)*I,dir(3)*I,scaleLoad,'LineWidth',2,'Color','b');
    1;
end


1;
end