function OUT = solve(model)
    % Wrapper. Calls a function stored in model.solver using the model as
    % argument
    OUT = model.solver(model);


end