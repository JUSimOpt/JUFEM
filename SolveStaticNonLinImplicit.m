function out = SolveStaticNonLinImplicit(model,varargin)
% Parameters:
% 'IterationConvergenceStudy' - 'On' or 'off'. Default 'Off'.   

%% Pre processing

IP = inputParser();
addParameter(IP,'IterationConvergenceStudy','off')
IP.parse(varargin{:});
IterationConvergenceStudy = strcmpi(IP.Results.IterationConvergenceStudy,'on');

% Solution field is a displacement vector
u0 = model.u0;

mesh = model.mesh;
steps = model.step;
dofs = model.dofs;

nodes = mesh.nodes;
points = mesh.P;

nnod = mesh.nnod; 
[nele,knod] = size(nodes);
neq = dofs*nnod; %Total number of equations
idofs = dofs*knod;

presc = model.presc; %Prescribed dofs
free = model.free; % Free dofs

% Choose ElementData function
UserElement = model.UserElement;
UserLoad = model.UserLoad;

[GP,GW] = mesh.IntegrationScheme(model.integrationOrder);

eqsInds = zeros(nele,idofs,'uint64');
switch dofs %Dimension specific parameters
    case 1
        error('1D code is not implemented!')
    case 2
        error('2D code is not implemented!')
    case 3
        for iel = 1:nele
            inods = nodes(iel,:);
            eqsInds(iel,1:3:end) = inods*3-2;
            eqsInds(iel,2:3:end) = inods*3-1;
            eqsInds(iel,3:3:end) = inods*3-0;
        end
end



materialData = model.materialData;

%% Pre-allocating output data
nSteps = length(steps);
maxnTimeIncs = 0;
elehistvar.ip(size(GP,1)).S = zeros(3); %Stress
elehistvar.ip(size(GP,1)).EE = zeros(3); %Elastic Strain
nodhistvar.u = zeros(neq,1);
nodhistvar.f = zeros(neq,1);
for istep = 1:nSteps
   step = steps(istep);
   analysisType = step.type;
   dTime = analysisType.initialTimeIncrement;
   totTime = analysisType.timePeriod;
   nTimeIncrements = totTime/dTime;
   for iinc = 1:nTimeIncrements
       out.step(nSteps).inc(iinc).elehistvar(nele) = elehistvar;
       out.step(nSteps).inc(iinc).nodhistvar = nodhistvar;
   end
   if nTimeIncrements > maxnTimeIncs
       maxnTimeIncs = nTimeIncrements;
   end
end


%% Main Solver
xfigure
hp = patch('Faces',mesh.Faces,'Vertices',mesh.P,'CData',zeros(nnod,1),'FaceColor','interp');
axis equal; view(3)

if IterationConvergenceStudy
    warning('Study must be run twice in order to get iteration errors.')
end
displog(model,'\n\n-------Static Implicit Non-Linear Solver-------\n');
displog(model,[TimeStamp(),'   Starting step\n']);tic1 = tic;
displog(model,'-----------------------------------------------\n');
txt0 =  sprintf('|  i  | Step |  time  |    dt   |  iter   |   r/r0     |\n');
txt1 =  sprintf('|-----|------|--------|---------|---------|------------|\n');
displog(model,txt0);
displog(model,txt1);
          
for istep = 1:nSteps
   % Loop through steps
   step = steps(istep);
   analysisType = step.type;
   maxIter = analysisType.maxEquilibriumIterations;
   dTime = analysisType.initialTimeIncrement;
   totTime = analysisType.timePeriod;
   nTimeIncrements = totTime/dTime;
   tol = analysisType.equilibriumTolerance;
   minTimeInc = analysisType.minTimeIncrement;
   
   if step.conservativeLoading==1
      % Compute the load vector here since we have conservative loading enabled
      % and thus the load vectors are independent of the deformed mesh.
      f = UserLoad(model,step);
   end
   
   
   u = model.u0; %Initialize the displacement field
   
   ctotIter = 0; %Total iterations
   cInc = 0;
   time = 0;
   while 1 % Loop over all time increments
      u_prev = u;
      time = time + dTime;
      cInc = cInc +1;
      if time > totTime
          break
      end
      
      loadFrac = time/totTime;
      f_frac = f*loadFrac;
      
      if IterationConvergenceStudy
          iterations = []; conditions = [];          
          try 
              u_exact = load('u_exact');
              u_exact = u_exact.u;
              if length(u_exact) ~= neq
                  u_exact = zeros(neq,1);
              end
              xfigure
              errs = [];hp_err = plot(0,0,'b-o'); hold on
              set(gca,'YScale','log')
              grid minor
              xlabel('Iteration')
              ylabel('\epsilon - displacement error')
              title(['Iteration convergence'])
              
          catch
              u_exact = zeros(neq,1);
          end
          xfigure
          hp1 = plot(0,0,'b-o'); hold on
          hp2 = plot(0,tol,'k-');
          xlabel('Iteration')
          ylabel('|R|/|R0| - Residual force')
          title(['Load fraction: ',num2str(loadFrac)])
      end
      
      if cInc > 1
          elehistvar_old =  out.step(istep).inc(cInc-1).elehistvar;
      else
          elehistvar_old = out.step(nSteps).inc(1).elehistvar;
      end
      elehistvar_new =  out.step(nSteps).inc(1).elehistvar;
      
      %Newton Iterations. Equilibrium Iterations
      for iter = 1:maxIter
          % Loop over all equilibrium iterations
          % Pre-allocate for assembler, enables paralization of the element
          % loop.    
          ctotIter = ctotIter + 1;
          
          
          g = zeros(neq,1); %Internal forces
          row = zeros(nele, idofs^2, 'double'); col = row; 
          valK = zeros(nele, idofs^2);
%           valr = zeros(nele, idofs);
          valg = zeros(nele, idofs);
          
          switch dofs %Dimension specific parameters
              case 1
                  error('1D code is not implemented!')
              case 2
                  error('2D code is not implemented!')
              case 3
                  % 3D
                  % Compute the loads here if we don't have conservative
                  % loading enabled. This means that the directions of the
                  % loads are dependent on the current geometry and need
                  % to be calculated. The Loading can also be a function of
                  % the continuum point.
                  if step.conservativeLoading==0
                     error('Non-linear loading not implemented. TODO') 
                  end
                  
                  
                  
                  % Computes the stiffness matrices, internal forces and
                  % resultant forces.
%                   disp('Assembling...'); tic;
                  parfor iel = 1:nele
                      % Loop over all elements
                      inods = nodes(iel,:);
                      ieqs = eqsInds(iel,:)';
                      ue = u(ieqs); %element displacements. OK overhead
%                       fe = f(ieqs);
                      Xc = points(inods,:);
                      
                      ielhistvar_old = elehistvar_old(iel);
                      
                      [Ke,ge,ielhistvar_new] = UserElement(ielhistvar_old,materialData,Xc,ue,ieqs,GP,GW,ctotIter,time,dTime,iter);
                      elehistvar_new(iel) = ielhistvar_new;
                      
                      % Assemble
                      ii = repmat(ieqs,1,idofs); jj = ii';
                      row(iel,:) = ii(:);
                      col(iel,:) = jj(:);
                      valK(iel,:) = Ke(:);
                      valg(iel,:) = ge(:);
                  end
                  K = sparse(row(:), col(:), valK(:), neq, neq);
                  % Assemble the resultant force vector
                  for iel = 1:nele
                      ieqs = eqsInds(iel,:)';
                      g(ieqs) = g(ieqs) + valg(iel,:)';
                  end
                  r = f_frac-g;
%                   toc
                  
%                   disp('Solving system...');tic;
                  % Solve the system
                  delta_u = u0;
                  delta_u(free) = K(free,free)\r(free);
%                   delta_u(free) = pcg(K(free,free),r(free));
%                   toc
                  % Update displacements
                  u = u + delta_u;
                  
                  out.step(istep).inc(cInc).elehistvar = elehistvar_new;
                  out.step(istep).inc(cInc).nodhistvar.u = u;
                  out.step(istep).inc(cInc).nodhistvar.f_frac = f_frac;
                  out.step(istep).inc(cInc).nodhistvar.g = g;
                  
                  
                  
                  % Iteration criteria
                  if(iter==1)
                      r0=r;
                  end
                  stopCond=norm(r(free))/norm(r0(free));
                  if stopCond > 1e5
%                      error('Iterations are diverging!') 
                    time = time - dTime;
                    cInc = cInc - 1;
                    dTime = dTime/2;
                    if dTime < minTimeInc
                       error('Min time increment reached!') 
                    end
                    u = u_prev;
                    break;
                  end
                  
                  
                  txt3 = sprintf('| %3d | %3d  | %0.4f | %0.4f  |  %3d    | %0.4e |  \n',ctotIter,istep,time,dTime,iter,stopCond);
                  displog(model,txt3);
                  
                  if IterationConvergenceStudy
                      iterations(end+1) = iter;
                      conditions(end+1) = stopCond;
                      set(hp1,'XData',iterations,'YData',conditions)
                      set(hp2, 'XData',iterations,'YData',conditions*0+tol)
                      
                      if ~all(u_exact==0)
                          errs(end+1) = norm(u-u_exact);
                          set(hp_err,'XData',iterations,'YData',errs)
                      end
                      
                      drawnow
                  end
                  
                  %Stop when resultant is less than tol of its original value
                  if (stopCond < tol)
                      break;
                  end
                  
                  %% Viz %TODO: Remove
                  U = [u(1:3:end), u(2:3:end), u(3:3:end)];
                  Ures = sum(U.^2,2).^0.5;
                  set(hp,'Vertices',mesh.P+U,'CData',Ures)
                  drawnow
          end
          
          if iter == maxIter
              error('Maximum equilibrium iterations exceeded. Probable divergence.')
          end
      end
      
%       if cInc > step.maxIncrement
%           error('Max increment reached.') %TODO: Figure out what Abaqus does here
%       end
      
      if IterationConvergenceStudy
          try
              save('u_exact','u')
          catch
              error('Unable to write "u_exact"')
          end
          if ~isempty(errs)
              ei1 = abs(log(errs(2:end)));
              ei0 = abs(log(errs(1:end-1)));
              rate = ei1./ei0;
              rate(end)=[];
              rate = rate'
              meanRate = mean(rate)
          end
      end
      
   end
   displog(model,[TimeStamp(),'\nStep done\n']);
   toc1 = toc(tic1);
   displog(model,['Elapsed time is ',ElapsedTime(toc1),'\n']);
    
end

out.model = model;

displog(model,[TimeStamp(),'\nAll complete!\n']);


end


function displog(model,msg)
    fprintf(msg)
    model.logfile.write(msg);
end

function str = TimeStamp()
    c = clock;
    str = sprintf('%d-%d-%d %d:%d:%2.1f',c(1),c(2),c(3),c(4),c(5),c(6));
end

function et = ElapsedTime(t)
    %t is in seconds
    nd = 0; nh = 0; nm = 0;
    et = '';
    if t >= 60*60*24 % more than a day
        nd = floor(t/(60*60*24));
        et = [num2str(nd),'d, '];
    end
    t = t - 60*60*24*nd;
    
    if t >= 60*60 % more than an hour
        nh = floor(t/3600);
        et = [et,num2str(nh),'h, '];
    end
    t = t - 60*60*nh;
    
    if t >= 60 % more than a minute
        nm = floor(t/60);
        et = [et,num2str(nm),'m, '];
    end
    t = t - 60*nm;
    
    et = [et,num2str(t),'s'];
    
end

