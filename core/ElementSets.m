classdef ElementSets
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        Type
        Nodes
    end
    
    methods
        function obj = ElementSets(name,type,nodes)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            % type is of classtype ElementType
            % TODO: Input validation
            obj.Name = name;
            obj.Type = type;
            obj.Nodes = nodes;
        end
        
    end
end

