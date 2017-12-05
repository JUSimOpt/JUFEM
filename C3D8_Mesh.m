classdef C3D8_Mesh
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    %     8-----7
    %    /|    /|
    %   5-----6 |
    %   | 4...|.3
    %   |/    |/ 
    %   1-----2
    % Faces are numbered in the following way
    % face1 = 1,2,3,4;
    % face1 = 5,6,7,8;
    % face1 = 1,2,6,5;
    % face1 = 2,3,7,6;
    % face1 = 3,4,8,7;
    % face1 = 4,1,5,8;
    
    
    
    properties
        nodes
        P
        NodeSets % contains Node sets
        ElementSets
        BoundarieConditions
        Loads 
        nele
        nnod
        Faces
    end
    
    properties (Hidden)
       XC
       YC
       ZC
    end
    
    methods
        function obj = C3D8_Mesh(P,nodes,varargin)
            %C3D8_Mesh Construct an instance of this class
            %   mesh = C3D8_Mesh(P,nodes)
            IP = inputParser;
            parse(IP,varargin{:});
            
            obj.P = P;
            obj.nodes = nodes;
            nele = size(nodes,1);
            obj.nele = nele;
            
            nnod = size(P,1);
            obj.nnod = nnod;
            
            faces = zeros(nele*6,4,'uint64');
            iface = [1,2,3,4;
                     5,6,7,8;
                     1,2,6,5;
                     2,3,7,6;
                     3,4,8,7;
                     4,1,5,8];
            lof = 1;
            for iel = 1:nele
                upf = lof+5;
                iv = nodes(iel,:);
                faces(lof:upf,:) = iv(iface);
                lof = upf +1;
            end
            obj.Faces = faces;
            
            obj.XC = P(:,1);
            obj.YC = P(:,2);
            obj.ZC = P(:,3);
            
            
        end
        
        function h = vizMesh(obj,varargin)
            %vizMesh Visualize the mesh in a figure
            %
            %   h = vizMesh()
            %   h = vizMesh(ele,properties)
            %   ele is a list of elements to display
            %   h contains visualization properties
            %            %
            %   properties:
            %   'NodeNumbers'       -   Display node numbers
            %   'ElementNumbers'    -   Display element numbers
            %
            
            ele = 1:size(obj.nodes,1);
            ele = ele(:);
            fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
            
            h.patch = patch(obj.XC(obj.Faces(fele(:),:)'),obj.YC(obj.Faces(fele(:),:)'),obj.ZC(obj.Faces(fele(:),:)'),'w','FaceColor','none');
            xlabel('X'); ylabel('Y'); zlabel('Z')
            axis equal tight
            set(gcf,'name','C3D8 Mesh')
            view(-55,45)
            
            if isenabled('NodeNumbers',varargin)
                if obj.nele > 100
                    warning('Cannot draw NodeNumbers, too many elements')
                    return
                end
                h.NodeText = [];
                %                 unique(T.Connectivity(:),'stable')'
                %                 1:T.nnod
                for i = 1:obj.nnod
                    h.NodeText = [h.NodeText; text(obj.XC(i),obj.YC(i),obj.ZC(i),num2str(i),'BackgroundColor','w') ];
                end
            end
            
            if isenabled('ElementNumbers',varargin)
                if obj.nele > 100
                    warning('Cannot draw ElementNumbers, too many elements')
                    return
                end
                
                h.EleText = [];
                for i = sort(ele)'
                    xm = mean(obj.XC(obj.Connectivity(i,:)));
                    ym = mean(obj.YC(obj.Connectivity(i,:)));
                    zm = mean(obj.ZC(obj.Connectivity(i,:)));
                    h.EleText = [h.EleText; text(xm,ym,zm,num2str(i),'BackgroundColor','y')];
                end
                
            end
            
            
        end










        
    end
 
end

