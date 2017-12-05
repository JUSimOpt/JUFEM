classdef AnalysisStep
    %AnalysisStep Step properties
    %   Used to set the step properties
    %
    % Properties:
    % 
    % name - String identifier
    % type - Type of the analysis, see StaticGeneralAnalyis class
    % maxIncrement - maximum number of increments in a step. Default 100
    % nonLinearGeometry - Default 0
    % solver - 'DirectSparse' (Default) or 'Iterative'
    
    properties
        name
        type
        maxIncrement
        nonLinearGeometry
        solver
        
    end
    
    methods
        function obj = AnalysisStep(name,type,varargin)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            IP = inputParser;
            addParameter(IP,'maxIncrement',100)
            addParameter(IP,'nonLinearGeometry',0)
            addParameter(IP,'solver','DirectSparse')
            parse(IP,varargin{:})
            obj.maxIncrement = IP.Results.maxIncrement;
            obj.nonLinearGeometry = IP.Results.nonLinearGeometry;
            obj.solver = IP.Results.solver;
            obj.name = name;
            obj.type = type;
        end
    end
end

