classdef NodeSets
    %NodeSets Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        Nodes
    end
    
    methods
        function obj = NodeSets(name,nodes)
            %NodeSets Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = name;
            obj.Nodes = nodes;
        end
    end
end

