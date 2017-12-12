classdef T3D2_Mesh
    %T3D2_Mesh Truss element in 3D with 2 nodes
    %   Detailed explanation goes here
    %
    %       ^n
    % ->    | 
    % o-----|-----o
    % 1           2  
    %

    
    
    
    properties
        nodes
        P
        NodeSets % contains Node sets
        ElementSets
        BoundarieConditions
        Loads 
        nele
        nnod
    end
    
    properties (Hidden)
       XC
       YC
       ZC
    end
    
    methods
        function obj = T3D2_Mesh(P,nodes,varargin)
            %T3D2_Mesh Construct an instance of this class
            %   mesh = T3D2_Mesh(P,nodes)
            IP = inputParser;
            parse(IP,varargin{:});
            
            obj.P = P;
            obj.nodes = nodes;
            nele = size(nodes,1);
            obj.nele = nele;
            
            nnod = size(P,1);
            obj.nnod = nnod;
                       
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
            
            error('Not Implemented! TODO')
            
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

        function [fi, B, detJ] = BaseFcnParam(obj,iel,iXi)
            % [fi, B, detJ] = BaseFcnParam(iel,iXi);
            %
            %       ^n
            % ->    | 
            % o-----|-----o
            % 1           2  
            %
            %
            % xi = iXi(1,:);
            %
            % xi in [0,1]
            locnods = obj.nodes(iel,:);
            Xc = obj.P(locnods,:);
            [fi, B, detJ] = Priv_BaseFcnParam(Xc,iXi);
        end
        
    end
    
    methods (Static)
        function [fi, B, detJ] = BaseFcnParam_Static(Xc,iXi)
            % [fi, B, detJ] = BaseFcnParam_Static(Xc,iXi)
            %
            %       ^n
            % ->    | 
            % o-----|-----o
            % 1           2  
            %
            % x = XC(:,1);
            % y = XC(:,2);
            % z = XC(:,3);
            %
            % xi = iXi(1,:);
            %
            % xi in [0,1]
            [fi, B, detJ] = Priv_BaseFcnParam(Xc,iXi);
        end
        
        function [fi, detJ] = BaseFcnParam_Static2(Xc,iXi)
            % [fi, detJ] = BaseFcnParam_Static(Xc,iXi)
            %
            %       ^n
            % ->    | 
            % o-----|-----o
            % 1           2  
            %
            % x = XC(:,1);
            % y = XC(:,2);
            % z = XC(:,3);
            %
            % xi = iXi(1,:);
            %
            % xi in [0,1]
            [fi, detJ] = Priv_BaseFcnParam2(Xc,iXi);
        end
        
        function [GP,GW] = IntegrationScheme(order)
            % [GP,GW] = IntegrationScheme(order)
            % Gauss Points and Weights on the domain [0,1]
            
            [GP,GW] = gauss(order, 0,1);
            GP = GP(:);
            GW = GW(:);
        end
        
    end
 
end

function [fi, B, detJ] = Priv_BaseFcnParam(XC,iXi)

xi = iXi(1);

X0 = XC(:,1);
Y0 = XC(:,2);
Z0 = XC(:,3);

fi = [1-xi, xi];
dfidxi = [-1,1];

t = XC(2,:)-XC(1,:);
n = [-t(2),t(1)];
n = n/norm(n);
detJ = norm(XC(2,:)-XC(1,:));

dxdxi = dfidxi*X0;
dydxi = dfidxi*Y0;
dzdxi = dfidxi*Z0;

if abs(dxdxi) < 1e-16
    dfidx = [0,0];
else
    dfidx = 1/dxdxi*dfidxi;
end

if abs(dydxi) < 1e-16
    dfidy = [0,0];
else
    dfidy = 1/dydxi*dfidxi;
end

if abs(dzdxi) < 1e-16
    dfidz = [0,0];
else
    dfidz = 1/dzdxi*dfidxi;
end

B = [dfidx;dfidy;dfidz];


end

function [fi, detJ] = Priv_BaseFcnParam2(XC,iXi)
xi = iXi(1);
fi = [1-xi, xi];
detJ = norm(XC(2,:)-XC(1,:));
end




