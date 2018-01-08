classdef StaticGeneralAnalyis
    %StaticGeneralAnalyis(parameter1,val1,parameter2,val2,...)
    % Parameters:
    %
    % initialTimeIncrement:
    % Initial time increment. This value will be modified as required if the 
    % automatic time stepping scheme is used or will be used as the constant 
    % time increment if the DIRECT parameter is used. If this entry is zero 
    % or is not specified, a default value that is equal to the total time 
    % period of the step is assumed. Default: 1
    %
    % timePeriod: Time period of the step. If this entry is zero or is not 
    % specified, a default value of  is assumed. Default: 1
    % 
    % minTimeIncrement: 
    % Minimum time increment allowed. Only used for automatic 
    % time incrementation. If Abaqus/Standard finds it needs a smaller time 
    % increment than this value, the analysis is terminated. If this entry is 
    % zero, a default value of the smaller of the suggested initial time 
    % increment or 1e-5 times the total time period is assumed.
    %
    % maxTimeIncrement:
    % Maximum time increment allowed. Only used for automatic time incrementation. 
    % If this value is not specified, no upper limit is imposed. Default: 1
    %
    % maxEquilibriumIterations:
    % upper limit on the number of consecutive equilibrium 
    % iterations. Default: 30
    %
    % equilibriumTolerance:
    % Used to stop the iterations when the resultant is less than equilibriumTolerance of its 
    % original value, i.e., |r|/|r0| < equilibriumTolerance
    %
    
    properties
        initialTimeIncrement
        timePeriod
        minTimeIncrement
        maxTimeIncrement
        maxEquilibriumIterations 
        equilibriumTolerance
    end
    
    methods
        function obj = StaticGeneralAnalyis(varargin)
            %StaticGeneralAnalyis()
            %StaticGeneralAnalyis(initialTimeIncrement, timePeriod, minTimeIncrement, maxTimeIncrement)
            IP = inputParser;
            addParameter(IP,'initialTimeIncrement',1)
            addParameter(IP,'timePeriod',1)
            addParameter(IP,'minTimeIncrement',1e-5)
            addParameter(IP,'maxTimeIncrement',1)
            addParameter(IP,'maxEquilibriumIterations',30)
            addParameter(IP,'equilibriumTolerance',1e-3)
            parse(IP,varargin{:})
            obj.initialTimeIncrement = IP.Results.initialTimeIncrement;
            obj.timePeriod = IP.Results.timePeriod;
            obj.minTimeIncrement = IP.Results.minTimeIncrement;
            obj.maxTimeIncrement = IP.Results.maxTimeIncrement;
            obj.maxEquilibriumIterations = IP.Results.maxEquilibriumIterations;
            obj.equilibriumTolerance = IP.Results.equilibriumTolerance;
        end
    end
end

