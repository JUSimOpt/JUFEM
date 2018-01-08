classdef M3D4_Mesh
    %M3D4_Mesh 4-node 3D Membrane mesh
    %   Detailed explanation goes here
    % 
    % 
    %   4-----3 
    %   |     |
    %   |     |
    %   1-----2

    
    
    
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
        function obj = M3D4_Mesh(P,nodes,varargin)
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
            if nargin > 1
               ele = varargin{2}; 
            end
            
            h.patch = patch('Faces',obj.nodes(ele,:),'Vertices',obj.P,'w','FaceColor','none');
            xlabel('X'); ylabel('Y'); zlabel('Z')
            axis equal tight
            set(gcf,'name','M3D4 Mesh')
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

        function [fi, detJ, B] = BaseFcnParam(obj,iel,iXi)
            % [fi, detJ, B] = BaseFcnParam(iel,iXi);
            % Unit cube with node one in 0,0
            %
            %   4-----3
            %   |     |
            %   |     |
            %   1-----2
            %
            %
            % xi = iXi(1,:);
            % eta = iXi(2,:);
            %
            % xi in [0,1]
            % eta in [0,1]
            locnods = obj.nodes(iel,:);
            Xc = obj.P(locnods,:);
            [fi, detJ, B] = Priv_BaseFcnParam(Xc,iXi);
        end
        
    end
    
    methods (Static)
        function [fi, detJ, B] = BaseFcnParam_Static(Xc,iXi)
            % [fi, detJ, B] = BaseFcnParam_Static(Xc,iXi)
            % Unit cube with node one in 0,0
            %
            %   4-----3
            %   |     |
            %   |     |
            %   1-----2
            %
            % x = XC(:,1);
            % y = XC(:,2);
            % z = XC(:,3);
            %
            % xi = iXi(1,:);
            % eta = iXi(2,:);
            %
            % xi in [0,1]
            % eta in [0,1]
            
            [fi, detJ, B] = Priv_BaseFcnParam(Xc,iXi);
        end
        
        function [GP,GW] = IntegrationScheme(order)
            % [GP,GW] = IntegrationScheme(order)
            % Gauss Points and Weights on the domain [0,1]
            
            [x,w] = gauss(order, 0,1);
            %x in [0,1]
            [X,Y]=meshgrid(x,x);
            GP = [X(:),Y(:)];
            
            [WX,WY]=meshgrid(w,w);
            GW = [WX(:),WY(:)];
            
            GW = prod(GW,2);
        end
        
        function [fi, B, Bxi, Beta, Bzeta, detJ] = BaseFcnParam_1Point_Static(XC0)
            error('Not implemented')
            % [fi, B, Bxi, Beta, Bzeta, detJ] = BaseFcnParam_1Point_Static(XC0)
            % Unit cube with node one in 0,0
            %
            %   4-----3
            %   |     |
            %   |     |
            %   1-----2
            %
            % xi in [0,1]
            % eta in [0,1]
            
            fi = [0.1250    0.1250    0.1250    0.1250    0.1250    0.1250    0.1250    0.1250];
            dfidxi = [-0.2500    0.2500    0.2500   -0.2500   -0.2500    0.2500    0.2500   -0.2500];
            dfideta = [-0.2500   -0.2500    0.2500    0.2500   -0.2500   -0.2500    0.2500    0.2500];
            dfidzeta = [-0.2500   -0.2500   -0.2500   -0.2500    0.2500    0.2500    0.2500    0.2500];
            
            
            X0 = XC0(:,1);
            Y0 = XC0(:,2);
            Z0 = XC0(:,3);
            
            J = [dfidxi*X0(:), dfidxi*Y0(:), dfidxi*Z0(:);...
                dfideta*X0(:), dfideta*Y0(:), dfideta*Z0(:);...
                dfidzeta*X0(:), dfidzeta*Y0(:), dfidzeta*Z0(:)];
            
            detJ = det(J);
            
            Bh = [dfidxi;dfideta;dfidzeta];
            
            Bhxi = [0,    0,    0,    0,    0,    0,   0,    0;
                0.5, -0.5,  0.5, -0.5,  0.5, -0.5, 0.5, -0.5;
                0.5, -0.5, -0.5,  0.5, -0.5,  0.5, 0.5, -0.5];
            
            Bheta = [0.5, -0.5,  0.5, -0.5,  0.5, -0.5, 0.5, -0.5;
                0,    0,    0,    0,    0,    0,   0,    0;
                0.5,  0.5, -0.5, -0.5, -0.5, -0.5, 0.5,  0.5];
            
            Bhzeta = [0.5, -0.5, -0.5,  0.5, -0.5,  0.5, 0.5, -0.5;
                0.5,  0.5, -0.5, -0.5, -0.5, -0.5, 0.5,  0.5;
                0,    0,    0,    0,    0,    0,   0,    0];
            
            BB = J\[Bh,Bhxi,Bheta,Bhzeta];
            B = BB(:,1:8);
            Bxi = BB(:,9:16);
            Beta = BB(:,17:24);
            Bzeta = BB(:,25:32);
            
            
        end
        
    end
 
end

function [fi, detJ, B] = Priv_BaseFcnParam(XC,iXi)

xi = iXi(1);
eta = iXi(2);

fi = [(1-xi)*(1-eta);
      xi*(1-eta);
      xi*eta;
      eta*(1-xi)].';


dfidxi = [eta - 1, 1 - eta, eta, -eta];

dfideta = [xi - 1, -xi, xi, 1 - xi];

dXdxi = dfidxi*XC;
dXdeta = dfideta*XC;


n = cross(dXdxi,dXdeta);
n = n/norm(n);

J = [dXdxi;...
    dXdeta;...
    n];

detJ = det(J);

B = J\[dfidxi;dfideta;0,0,0,0];

end


function [x, w, A] = gauss(n, a, b)

%------------------------------------------------------------------------------
% gauss.m
%------------------------------------------------------------------------------
%
% Purpose:
%
% Generates abscissas and weigths on I = [ a, b ] (for Gaussian quadrature).
%
%
% Syntax:
%
% [x, w, A] = gauss(n, a, b);
%
%
% Input:
%
% n    integer    Number of quadrature points.
% a    real       Left endpoint of interval.
% b    real       Right endpoint of interval.
%
%
% Output:
%
% x    real       Quadrature points.
% w    real       Weigths.
% A    real       Area of domain.
%------------------------------------------------------------------------------


% 3-term recurrence coefficients:
n = 1:(n - 1);
beta = 1 ./ sqrt(4 - 1 ./ (n .* n));

% Jacobi matrix:
J = diag(beta, 1) + diag(beta, -1); 


% Eigenvalue decomposition:

%
% e-values are used for abscissas, whereas e-vectors determine weights.
%

[V, D] = eig(J);
x = diag(D);


% Size of domain:
A = b - a;


% Change of interval:

%
% Initally we have I0 = [ -1, 1 ]; now one gets I = [ a, b ] instead.
%
% The abscissas are Legendre points.
%

if ~(a == -1 && b == 1)
  x = 0.5 * A * x + 0.5 * (b + a);
end


% Weigths:
w = V(1, :) .* V(1, :);
end
