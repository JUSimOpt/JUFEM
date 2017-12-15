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
numberOfRuns = 1; 

mesh = CreateMesh(modelFileName, refinements);
mesh = CreateSets(mesh);
mesh = CreateBoundaryConditions(mesh);
mesh = CreateLoads(mesh);

% debug figure
xfigure; hold on;
h = mesh.vizMesh;
h.patch.FaceColor = 'w';
light    
view(20, -320)

%% creating locking plot
% finding similar node as used in Peters paper. 
solNode = find(roundNR(mesh.P(:,1),4) == roundNR(max(mesh.P(:,1)),4) ...
    & roundNR(mesh.P(:,2),4) == roundNR(min(mesh.P(:,2)),4) ...
    & roundNR(mesh.P(:,3),4) == roundNR(max(mesh.P(:,3)),4));
solIeqs = 3*solNode;

scatter3(mesh.P(solNode,1), mesh.P(solNode,2), mesh.P(solNode,3),100,'r','filled')
hold off

nu_StoreValues = zeros(numberOfRuns,1);
u_StoreValues = nu_StoreValues; 
meanRate = nu_StoreValues; 
nu_Vector = linspace(0.45,0.4999,numberOfRuns)';
for iRuns = 1:numberOfRuns
%     close all
%     clear model    
    % materialData = CreateMaterialData();
    % Material properties
    E = 210000; %MPa
    
    nu = nu_Vector(iRuns);
    nu = 0;
    K=E/3/(1-2*nu);
    mu=E/2/(1+nu);
    mu1 = 0.5*mu;
    mu2 = mu-mu1;
    materialData(1).Type = 'Neo-Hooke';
    materialData(1).materialFcn = @(gradU)NeoHook3D_2PK(gradU,mu1/2,K);
    model.dofs = 3;
    procedureType = StaticGeneralAnalyis(   'initialTimeIncrement', 1, ...
        'timePeriod', 1, ...
        'minTimeIncrement', 1/100, ...
        'maxTimeIncrement', 1,...
        'equilibriumTolerance', 1e-3,...
        'maxEquilibriumIterations', 30);
    
    model.step(1) = AnalysisStep('Step-1', procedureType);
    model.step(1).nonLinearGeometry = 1;
    model.step(1).conservativeLoading = 1;
    
    % model.UserElement = @ElementDataFullIntegration_3D; %Default
    model.UserElement = @ElementData1PointIntegration_3D; %Default
    model.integrationOrder = 2;
    
    model.UserLoad = @StaticLoad; %Default
    model.solver = @(model)SolveStaticNonLinImplicit(model,'IterationConvergenceStudy','on');
    % model.baseFcnParam = @mesh.BaseFcnParam_Static; %Default
    model.baseFcnParam = @mesh.BaseFcnParam_1Point_Static; % One point integration
    
    %% Pre-processing model
    % Runs an input validation and optionally writes and input file
    % TODO inport model from .mat and run through the PreProcessor
    model = PreProcessor(mesh,model,materialData);
    
    vizBCandLoads(model);
    
    %% Solver be here
    OUT = model.solver(model);
    nu_StoreValues(iRuns) = nu; 
    u_StoreValues(iRuns) = OUT.u(solIeqs);
    meanRate(iRuns) = OUT.meanRate(iRuns);
end
save('nu_vals', 'nu_StoreValues')
save('u_1pointIntegration', ' u_StoreValues')
% save('u_fullIntegration', ' u_StoreValues')

%% creating Locking plot
nuVals = load('nu_vals');
uOnePoint = load('u_1pointIntegration');
% uFullInteg = load('u_fullIntegration');
xfigure; hold on; grid on; axis equal;
xlabel('Poisson Ratio, \nu')
ylabel('Displacement')
% plot(nuVals,uOnePoint,'-ko', 'LineWidth',1, 'MarkerSize',5, 'MarkerEdgeColor','k')
% plot(nuVals,uFullInteg,'--ko', 'LineWidth',1, 'MarkerSize',5, 'MarkerEdgeColor','k')

function mesh = CreateMesh(inFile, refinements) 
%% Run gmsh
% User defined  
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

% Node sets (user defined)
mesh.NodeSets(1).name = 'Set-BC';
mesh.NodeSets(1).nodes = find(mesh.P(:,1) == min(mesh.P(:,1)));

mesh.NodeSets(2).name = 'SurfaceLoadSet';
mesh.NodeSets(2).nodes = FindSurfaceNodes(mesh, max(mesh.P(:,1)));
end

function mesh = CreateBoundaryConditions(mesh)
% User defined
mesh.BoundarieConditions(1).Name = 'Clamp';
mesh.BoundarieConditions(1).Type = BCType.DispRot;
mesh.BoundarieConditions(1).NodeSet = mesh.NodeSets(1).nodes;
mesh.BoundarieConditions(1).U1 = 0;
mesh.BoundarieConditions(1).U2 = 0;
mesh.BoundarieConditions(1).U3 = 0;
end

function mesh = CreateLoads(mesh)
% User Defined
%TODO: Select by name somehow instead of indices
mesh.Loads(1).Name = 'TractionSurfaceLoad';
mesh.Loads(1).Type = LoadType.SurfaceLoad;
mesh.Loads(1).Magnitude = 1000; %N/mm
mesh.Loads(1).Direction = [0,0,-1];
mesh.Loads(1).Set = mesh.NodeSets(2).nodes;
mesh.Loads(1).IntegrationOrder = 2; %Linear load on edge
end




