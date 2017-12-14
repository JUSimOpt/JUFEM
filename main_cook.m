clear,clc,close all

% This is a demo file solving a Hyperelastic model
% Models are created using the user defined functions in this file. 
% The mesh is created using gmsh and an .geo file describing the geometry.
% TODO: Create an input file using ANSA
% 
%
% Add the necessary packages
addpath(fullfile('Gmsh','Gmsh'))
addpath(fullfile('Hex1Mesh','Hex1Mesh'))
addpath(fullfile('xFigure'))

%% Create Model
modelFileName = 'CooksMembrane.geo';

refinements = 1;

mesh = CreateMesh(modelFileName, refinements);
% debug figure
xfigure
h = mesh.vizMesh;
h.patch.FaceColor = 'w';
light    
view(20, -320)


materialData = CreateMaterialData();
model.dofs = 3;
procedureType = StaticGeneralAnalyis(   'initialTimeIncrement', 1/2, ...
                                        'timePeriod', 1, ...
                                        'minTimeIncrement', 1/100, ...
                                        'maxTimeIncrement', 1,...
                                        'equilibriumTolerance', 1e-3,...
                                        'maxEquilibriumIterations', 30);
                                    
model.step(1) = AnalysisStep('Step-1', procedureType);
model.step(1).nonLinearGeometry = 1;

mesh = CreateSets(mesh);
mesh = CreateBoundaryConditions(mesh);
mesh = CreateLoads(mesh);
model = PreProcessor(mesh,model,materialData);

vizBCandLoads(model);

%% Solver be here
OUT = solve(model);

%% User Defined functions for creating the model
function materialData = CreateMaterialData()

type = 'Neo-Hooke';
E = 210000; %Yungs Modulus MPa
nu = 0.3333; %Poissons Ratio

K = E/3/(1-2*nu); %Bulk Modulus
mu=E/2/(1+nu);


materialData(1).Type = type;
materialData(1).E = E;
materialData(1).nu = nu;
materialData(1).K = K;
materialData(1).mu = mu;
materialData(1).materialFcn = @NeoHook3D_2PK;

materialData(1).K = E/3/(1-2*nu); %Bulk Modulus
end

function mesh = CreateMesh(inFile, refinements)
% User defined
%     x0=0;x1=10;
%     y0=0;y1=2;
%     z0=0;z1=2;
%     ratio = round((x1-x0)/(y1-y0));
%     nxe = ratio*refinements;
%     nye = refinements;
%     nze = refinements;
%     mesh = HexP1MeshAbaqus(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
    

%% Run gmsh
% scanning current document
fid = fopen(inFile,'r');
if fid == -1
    error(['Cannot open the file: ', inFile])
end
s = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

% changing to current refinement
size_s = size(s{1},1);
s = s{1};
for refine = 1:refinements
    size_s = size_s+2;
    s{size_s} = '';
    s{size_s} = [s{size_s},'RefineMesh;'];
end
save('TempMatlab_refined.geo','s');
newInFile = 'TempMatlab_refined.geo.geo';

fid = fopen(newInFile,'w');
try
    fprintf(fid,'%s\n',s{:});
catch
    fclose(fid);
    error('Something went wrong writing to the file!');
end
fclose(fid);

% scanning document and generating mesh with current refinement
[status,cmdout] = RunGmshScript(newInFile,'verbose','on');

% extract element types 3 and 5 (Quads and Hex)
msh = MshRead(fullfile(pwd,'mesh.msh'),'typesToExtract',[3,5]);

% remap Gmsh hex elements to Abaqus numbering
map = [5,6,2,1,8,7,3,4];
mesh = C3D8_Mesh(msh.P,msh.ElementList(2).nodes(:,map));

end

function mesh = CreateSets(mesh)
% User defined
% TODO: Read from inp file.

%% Node sets (user defined)
mesh.NodeSets(1).name = 'Set-BC';
mesh.NodeSets(1).nodes = find(mesh.P(:,1) == min(mesh.P(:,1)));

% mesh.NodeSets(2).name = 'PointLoadSet';
% mesh.NodeSets(2).nodes = find(  mesh.P(:,1) == max(mesh.P(:,1)) & ...
%                                 mesh.P(:,2) == min(mesh.P(:,2)) & ... 
%                                 mesh.P(:,3) == max(mesh.P(:,3)) );

mesh.NodeSets(2).name = 'SurfaceLoadSet';
mesh.NodeSets(2).nodes  = findValInVector(mesh.P(:,1), max(mesh.P(:,1)));

%% Element sets
mesh.ElementSets(1).name = 'All-Hex';
mesh.ElementSets(1).type = ElementType.C3D8R;
mesh.ElementSets(1).nodes = mesh.nodes;

mesh.ElementSets(2).name = 'LoadElementSet-1';
mesh.ElementSets(2).type = ElementType.TwoNode3D_Line;

% Define line elements at x=x1 and z = z1
nod = find(mesh.P(:,1)==max(mesh.P(:,1)) & mesh.P(:,3)==max(mesh.P(:,3)));
ele = find(any(ismember(mesh.nodes,nod),2));
nodes = NaN(length(ele),2);
for i = 1:length(ele)
    iel = ele(i);
    iv = mesh.nodes(iel,:);
    nodes(i,:)=iv(ismember(iv,nod));
end
mesh.ElementSets(2).nodes = nodes;

mesh.ElementSets(3).name = 'LoadElementSet-2';
mesh.ElementSets(3).type = ElementType.TwoNode3D_Line;

% Define line elements at x=x1 and z = z0
nod = find(mesh.P(:,1)==max(mesh.P(:,1)) & mesh.P(:,3)==min(mesh.P(:,3)));
ele = find(any(ismember(mesh.nodes,nod),2));
nodes = NaN(length(ele),2);
for i = 1:length(ele)
    iel = ele(i);
    iv = mesh.nodes(iel,:);
    nodes(i,:)=iv(ismember(iv,nod));
end
mesh.ElementSets(3).nodes = nodes;



end

function mesh = CreateBoundaryConditions(mesh)
% User defined

mesh.BoundarieConditions(1).Name = 'Clamp';
mesh.BoundarieConditions(1).Type = BCType.DispRot;
mesh.BoundarieConditions(1).NodeSet = mesh.NodeSets(1).nodes;
mesh.BoundarieConditions(1).U1 = 0;
mesh.BoundarieConditions(1).U2 = 0;
mesh.BoundarieConditions(1).U3 = 0;
% mesh.BoundarieConditions(1).UR1 = 0;
% mesh.BoundarieConditions(1).UR2 = 0;
% mesh.BoundarieConditions(1).UR3 = 0;






end

function mesh = CreateLoads(mesh)
% User Defined

ntype = 2;
%TODO: Select by name somehow instead of indices
switch ntype
    case 1
        mesh.Loads(1).Name = 'UpperLineLoad';
        mesh.Loads(1).Type = LoadType.EdgeLoad;
        mesh.Loads(1).Magnitude = 0.4; %N/mm
        mesh.Loads(1).Direction = [-1,0,0];
        mesh.Loads(1).Set = mesh.ElementSets(2).nodes;
        
        mesh.Loads(2).Name = 'LowerLineLoad';
        mesh.Loads(2).Type = LoadType.EdgeLoad;
        mesh.Loads(2).Magnitude = 0.4; %N/mm
        mesh.Loads(2).Direction = [1,0,0];
        mesh.Loads(2).Set = mesh.ElementSets(3).nodes;
        
        mesh.Loads(3).Name = 'PointLoad';
        mesh.Loads(3).Type = LoadType.PointLoad;
        mesh.Loads(3).Magnitude = 0.2; %N
        mesh.Loads(3).Direction = [0,1,0];
        mesh.Loads(3).Set = mesh.NodeSets(2).nodes;
    case 2 % surface load        
        mesh.Loads(1).Name = 'TractionSurfaceLoad';
        mesh.Loads(1).Type = LoadType.SurfaceLoad;
        mesh.Loads(1).Magnitude = 1; %N/mm
        mesh.Loads(1).Direction = [0,-1,0];
        mesh.Loads(1).Set = mesh.NodeSets(2).nodes;
end

end

function dir = UpDir(dir,n)
    for i = 1:n
        dir = fileparts(dir);
    end
end



