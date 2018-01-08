clear,clc,close all

% This is a demo file solving a Hyperelastic model
% Models are created using the user defined functions in this file. 
% The mesh is created using gmsh and a .geo file describing the geometry.
% TODO: Create an input file using ANSA
% 
%
% Add the necessary packages
addpath(fullfile('core')) %All JUFEM functions are located here
addpath(fullfile('core','Gmsh','Gmsh'))
addpath(fullfile('core','xFigure'))





%% Create Model - User defined
model.name = 'Reese';

refinements = 0; d = 3; thin = 0;
intRules = {'Full','Reduced','IsoStab','IsoVolStab'};
intRule = intRules{1};

for d = 4:4
    for j = 1:4
        intRule = intRules{j}

mesh = CreateMesh(refinements,d,thin);
mesh = CreateSets(mesh);
mesh = CreateBoundaryConditions(mesh);
mesh = CreateLoads(mesh);

xfigure
h = mesh.vizMesh;
h.patch.FaceColor = 'w';
light
1;



% Material properties
E = 200; %GPa
nu = 0.3;
K=E/(3*(1-2*nu));
mu=E/(2*(1+nu));
mu1 = 0.5*mu;
mu2 = mu-mu1;
materialData(1).Type = 'Neo-Hooke';
materialData(1).materialFcn = @(gradU)NeoHook3D_2PK(gradU,mu1,K);

procedureType = StaticGeneralAnalyis('initialTimeIncrement',1/3,'timePeriod',1,'minTimeIncrement',1/100,'maxTimeIncrement',1,...
                                     'equilibriumTolerance',1e-3,'maxEquilibriumIterations',30);
model.step(1) = AnalysisStep('Step-1',procedureType);
model.step(1).nonLinearGeometry = 1;
model.step(1).conservativeLoading = 1;


switch intRule
    case 'Full'
model.UserElement = @(histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter)UserElement_3D_Full(...
                     histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter,...
                     @mesh.BaseFcnParam_Static ); % Default
    case 'Reduced'             
model.UserElement = @(histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter)UserElement_3D_Reduced(...
                     histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter,...
                     @mesh.BaseFcnParam_Static);
    case 'IsoStab'  
model.UserElement = @(histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter)UserElement_3D_1Point_isoStab(...
                     histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter,...
                     @mesh.BaseFcnParam_Static,@mesh.BaseFcnParam_1Point_Static );
    case 'IsoVolStab'
model.UserElement = @(histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter)UserElement_3D_1Point_IsoVolStab(...
                     histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter,...
                     @mesh.BaseFcnParam_Static,@mesh.BaseFcnParam_1Point_Static );  
end

model.integrationOrder = 2;

model.UserLoad = @StaticLoad; %Default
model.solver = @(model)SolveStaticNonLinImplicit(model,'IterationConvergenceStudy','off');



%% Pre-processing model
% Runs an input validation and optionally writes and input file
% TODO inport model from .mat and run through the PreProcessor
model = PreProcessor(mesh,model,materialData);

vizBCandLoads(model);

%% Solver be here
OUT = model.solver(model);

% filename = ['results//ReeseBeam_Thick_d',num2str(d),'_',intRule,'.mat']
% save(filename,'model','OUT')

    end
end

%% User Defined functions for creating the model

function mesh = CreateMesh(refinements,d,thin)
% User defined
    
    %% Modify .geo file
    inFile = 'ReeseBeam.geo';
    fid = fopen(inFile,'r');
    if fid == -1
        error(['Cannot open the file: ', inFile])
    end
    try
        s = textscan(fid,'%s','Delimiter','\n');
    catch
        fclose(fid);
        error('Something went wrong reading the file!');
    end
    fclose(fid);
    s = s{1};
    
    if thin
        s{7} = 'W=0.2; //mm';
        s{8} = 'H=0.2; //mm';
    else
        s{7} = 'W=2; //mm';
        s{8} = 'H=2; //mm';
    end
    
    s{9} = ['d = ',num2str(d),'; //mm'];
    
    RefineMeshLine = 119;
    s{RefineMeshLine} = '';
    for i = 1:refinements
        s{RefineMeshLine} = [s{RefineMeshLine},'RefineMesh;'];
    end
    
    fid = fopen(inFile,'w');
    if fid == -1
        error(['Cannot open the file: ', inFile])
    end
    try
        fprintf(fid,'%s\n',s{:});
    catch
        fclose(fid);
        error('Something went wrong writing to the file!');
        
    end
    fclose(fid);

    %% Run gmsh
    [status,cmdout] = RunGmshScript(inFile,'verbose','on','OutFile','mesh.msh');
    msh = MshRead(fullfile(pwd,'mesh.msh'),'typesToExtract',[3,5]); %Extract element types 3 and 5 (Quads and Hex)
    1;
    %remap Gmsh hex elements to Abaqus numbering
    map = [5,6,2,1,8,7,3,4];
    mesh = C3D8_Mesh(msh.P,msh.ElementList(2).nodes(:,map));
    
end

function mesh = CreateSets(mesh)
% User defined
% TODO: Read from inp file.

%% Node sets
mesh.NodeSets(1).name = 'Set-BC'; % User defined
mesh.NodeSets(1).nodes = find(mesh.P(:,1)==min(mesh.P(:,1)));

mesh.NodeSets(2).name = 'PointLoadSet';
mesh.NodeSets(2).nodes = find(mesh.P(:,1)==max(mesh.P(:,1)) & mesh.P(:,2)==min(mesh.P(:,2)) & mesh.P(:,3)==max(mesh.P(:,3)) );

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

%TODO: Select by name somehow instead of indices

mesh.Loads(1).Name = 'UpperLineLoad';
mesh.Loads(1).Type = LoadType.EdgeLoad;
% mesh.Loads(1).Magnitude = 0.4; %N/mm
mesh.Loads(1).Magnitude = 4/10; %N/mm
mesh.Loads(1).Direction = [-1,0,0];
mesh.Loads(1).Set = mesh.ElementSets(2).nodes;
mesh.Loads(1).IntegrationOrder = 1; %Linear load on edge

mesh.Loads(2).Name = 'LowerLineLoad';
mesh.Loads(2).Type = LoadType.EdgeLoad;
% mesh.Loads(2).Magnitude = 0.4; %N/mm
mesh.Loads(2).Magnitude = 4/10; %N/mm
mesh.Loads(2).Direction = [1,0,0];
mesh.Loads(2).Set = mesh.ElementSets(3).nodes;
mesh.Loads(2).IntegrationOrder = 1; %Linear load on edge

mesh.Loads(3).Name = 'PointLoad';
mesh.Loads(3).Type = LoadType.PointLoad;
% mesh.Loads(3).Magnitude = 0.2; %N
mesh.Loads(3).Magnitude = 0.03333; %N/mm
mesh.Loads(3).Direction = [0,1,0];
mesh.Loads(3).Set = mesh.NodeSets(2).nodes;


end





