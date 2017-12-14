classdef AnalysisStep
    %AnalysisStep Step properties
    %   Used to set the step properties
    %
    % Parameters:
    % 
    % name - String identifier
    % type - Type of the analysis, see StaticGeneralAnalyis class
    % maxIncrement - maximum number of increments in a step. Default 100
    % nonLinearGeometry - Default: 0
    % conservativeLoading - Default: 1
    % solver - 'DirectSparse' (Default) or 'Iterative'
    
    properties
        name
        type
        maxIncrement
        nonLinearGeometry
        conservativeLoading
        solver
        
    end
    
    methods
        function obj = AnalysisStep(name,type,varargin)
            %AnalysisStep Construct an instance of this class
            %   AnalysisStep(name,type,Parameters)
            % Parameters:
            % 
            % name - String identifier
            % type - Type of the analysis, see StaticGeneralAnalyis class
            % maxIncrement - maximum number of increments in a step. Default 100
            % nonLinearGeometry - Default: 0
            % convervativeLoading - Default: 1
            % solver - 'DirectSparse' (Default) or 'Iterative'
            
            IP = inputParser;
            addParameter(IP,'maxIncrement',100)
            addParameter(IP,'nonLinearGeometry',0)
            addParameter(IP,'solver','DirectSparse')
            addParameter(IP,'conservativeLoading',1)
            parse(IP,varargin{:})
            obj.maxIncrement = IP.Results.maxIncrement;
            obj.nonLinearGeometry = IP.Results.nonLinearGeometry;
            obj.solver = IP.Results.solver;
            obj.conservativeLoading = IP.Results.conservativeLoading;
            obj.name = name;
            obj.type = type;
        end
    end
end

