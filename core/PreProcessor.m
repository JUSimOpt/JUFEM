function model = PreProcessor(mesh,model,materialData)
% Processor
% Everything is contained in the model struct.
% TODO: Write to ascii file. .inp possibly in order to run in abaqus.

%% Input Validation
% TODO: Don't trust the user, check for logical constistancy. Throw error
% messages if needed.


model.dofs = size(mesh.P,2);

model = processBC(mesh,model);
model.materialData = materialData;


%% Write log file
filename = [model.name,'.log'];
model.logfile = LogFile(filename,'w');
s = writeLog(model);


%% Write inp file
% TODO:

end



















function outString = writeLog(model)
outString = [];
outString = [outString,sprintf('--------PreProcessor---------\n')];
outString = [outString,sprintf('------------------------------------\n\n')];
outString = [outString,sprintf('* Mesh info:\n')];
outString = [outString,sprintf('Number of elements: %d\n',model.mesh.nele)];
outString = [outString,sprintf('Number of nodes: %d\n',model.mesh.nnod)];
outString = [outString,sprintf('Degrees of freedom per node: %d\n',model.dofs)];
outString = [outString,sprintf('Total degrees of freedom: %d\n\n',length(model.u0))];
outString = [outString,sprintf('* Material info:\n')];
fnames = fieldnames(model.materialData);
for i = 1:length(fnames)
    val = model.materialData.(fnames{i});
    if isa(val,'function_handle')
        outString = [outString,sprintf('Material.%s: %s\n',fnames{i},func2str(val))];
    else
        outString = [outString,sprintf('Material.%s: %s\n',fnames{i},num2str(val))];
    end
    
end
outString = [outString,sprintf('\n* Model info:\n')];

nsteps = length(model.step);
for i = 1:nsteps
    outString = [outString,sprintf('Step: %s\n',model.step(i).name)];
    outString = [outString,sprintf('Analysis type: %s\n',class(model.step(i).type))];
    outString = [outString,sprintf('nonLinearGeometry: %d\n',model.step(i).nonLinearGeometry)];
    outString = [outString,sprintf('maxIncrement: %d\n',model.step(i).maxIncrement)];
    outString = [outString,sprintf('System solver: %s\n',model.step(i).solver)];
end

outString = [outString,sprintf('\n*Solver: %s\n',func2str(model.solver))];

fprintf(outString)


model.logfile.write(outString);
model.logfile.permission = 'a'; % Set logfile to append


end

function model = processBC(mesh,model)
dofs = model.dofs; % TODO: where should this go?
neq = mesh.nnod*dofs;
u0 = zeros(neq,1);

nb = length(mesh.BoundarieConditions);
presc=[];
for i = 1:nb
    ind=mesh.BoundarieConditions(i).NodeSet;
    presc_x = []; presc_y=presc_x; presc_z=presc_x;
    if isfield(mesh.BoundarieConditions(i),'U1')
        presc_x = 3*ind-2; % Prescribed x
        u0(presc_x) = 0;
    end
    if isfield(mesh.BoundarieConditions(i),'U2')
        presc_y = 3*ind-1; % Prescribed y
        u0(presc_y) = 0;
    end
    if isfield(mesh.BoundarieConditions(i),'U3')
        presc_z = 3*ind-0; % Prescribed z
        u0(presc_z) = 0;
    end
    
    mesh.BoundarieConditions(i).presc_x=presc_x;
    mesh.BoundarieConditions(i).presc_y=presc_y;
    mesh.BoundarieConditions(i).presc_z=presc_z;
    
    presc=[presc;presc_x;presc_y;presc_z]; % TODO: Fix this by preallocating
end
free = setdiff(1:neq,presc)';


model.u0 = u0;
model.presc = presc;
model.free = free;
model.mesh = mesh;

end