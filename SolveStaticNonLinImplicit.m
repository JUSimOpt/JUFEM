function OUT = SolveStaticNonLinImplicit(model)
1



steps = model.step;
for istep = 1:length(steps)
   % Lopp through steps
   step = steps(istep);
   analysisType = step.type;
   nTimeIncrements = analysisType.timePeriod/analysisType.initialTimeIncrement;
   cInc = 1; %Increment counter
   for iTime = 1:nTimeIncrements
      % Loop over all time increments
      
      
      for iter = 1:analysisType.maxEquilibriumIterations
          % Loop over all equilibrium iterations
          
          
          
          
          if iter == analysisType.maxEquilibriumIterations
              error('Maximum equilibrium iterations exceeded. Probable divergence.')
          end
      end
      
      if cInc > step.maxIncrement
          error('Max increment reached.') %TODO: Figure out what Abaqus doeas here
      end
      cInc = cInc + 1;
   end
    
end






% Newton Iteration

% Element Loop

% Material

% Assemble

% Solve linear system

OUT = 1;
end



