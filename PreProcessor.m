function model = PreProcessor(mesh,model,materialData)
% Processor
% Everything is contained in the model struct.
% TODO: Write to ascii file. .inp possibly in order to run in abaqus.

model = processBC(mesh,model);
model.materialData = materialData;

s = writeLog(model);

end

function outString = writeLog(model)
outString = [];
outString = [outString,sprintf('--------PreProcessor---------\n')];
outString = [outString,sprintf('------------------------------------\n\n')];
outString = [outString,sprintf('* Mesh info:\n')];
outString = [outString,sprintf('Number of elements: %d\n',model.mesh.nele)];
outString = [outString,sprintf('Number of nodes: %d\n',model.mesh.nnod)];
outString = [outString,sprintf('Degrees of freedom per node: %d\n',model.dofs)];
outString = [outString,sprintf('Total degrees of freedom: %d\n\n',length(model.u))];
outString = [outString,sprintf('* Material info:\n')];
fnames = fieldnames(model.materialData);
for i = 1:length(fnames)
    val = model.materialData.(fnames{i});
    if isa(val,'function_handle')
        outString = [outString,sprintf('Material.%s: %s\n',fnames{i},func2str(val))];
    else
        outString = [outString,sprintf('Material.%s: %s\n',fnames{i},val)];
    end
    
end
outString = [outString,sprintf('\n* Model info:\n')];
outString = [outString,sprintf('nIncrements: %d\n',model.nIncrements)];
outString = [outString,sprintf('maxiterations: %d\n',model.maxIterations)];
outString = [outString,sprintf('Newton tolerance: %s\n',model.tol)];
outString = [outString,sprintf('Solver: %s\n',func2str(model.solver))];

fprintf(outString)

fID = fopen('preProcessor.log','w');
if fID == -1
    error('Unable to write to file "preProcessor.log"')
end

try
    fprintf(fID,outString);
catch
    fclose(fID);
    error('Error writing to file "preProcessor.log"')
end
fclose(fID);


end

function model = processBC(mesh,model)
dofs = model.dofs; % TODO: where should this go?
neq = mesh.nnod*dofs;
u = zeros(neq,1);

nb = length(mesh.BoundarieConditions);
presc=[];
for i = 1:nb
    ind=mesh.BoundarieConditions(i).NodeSet;
    presc_x = []; presc_y=presc_x; presc_z=presc_x;
    if isfield(mesh.BoundarieConditions(i),'U1')
        presc_x = 3*ind-2; % Prescribed x
        u(presc_x) = 0;
    end
    if isfield(mesh.BoundarieConditions(i),'U2')
        presc_y = 3*ind-1; % Prescribed y
        u(presc_y) = 0;
    end
    if isfield(mesh.BoundarieConditions(i),'U3')
        presc_z = 3*ind-0; % Prescribed z
        u(presc_z) = 0;
    end
    
    mesh.BoundarieConditions(i).presc_x=presc_x;
    mesh.BoundarieConditions(i).presc_y=presc_y;
    mesh.BoundarieConditions(i).presc_z=presc_z;
    
    presc=[presc;presc_x;presc_y;presc_z]; % TODO: Fix this by preallocating
end
free = setdiff(1:neq,presc)';


model.u = u;
model.presc = presc;
model.free = free;
model.mesh = mesh;

end