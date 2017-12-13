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

        function [fi, detJ, B] = BaseFcnParam(obj,iel,iXi)
            % [fi, detJ, B] = BaseFcnParam(iel,iXi);
            % Unit cube with node one in 0,0,0
            %     8-----7
            %    /|    /|
            %   5-----6 |
            %   | 4...|.3
            %   |/    |/
            %   1-----2
            %
            %
            % xi = iXi(1,:);
            % eta = iXi(2,:);
            % zeta = iXi(3,:);
            %
            % xi in [0,1]
            % eta in [0,1]
            % zeta in [0,1]
            locnods = obj.nodes(iel,:);
            Xc = obj.P(locnods,:);
            [fi, detJ, B] = Priv_BaseFcnParam(Xc,iXi);
        end
        
    end
    
    methods (Static)
        function [fi, detJ, B] = BaseFcnParam_Static(Xc,iXi)
            % [fi, detJ, B] = BaseFcnParam_Static(Xc,iXi)
            % Unit cube with node one in 0,0,0
            %     8-----7
            %    /|    /|
            %   5-----6 |
            %   | 4...|.3
            %   |/    |/
            %   1-----2
            %
            % x = XC(:,1);
            % y = XC(:,2);
            % z = XC(:,3);
            %
            % xi = iXi(1,:);
            % eta = iXi(2,:);
            % zeta = iXi(3,:);
            %
            % xi in [0,1]
            % eta in [0,1]
            % zeta in [0,1]
            [fi, detJ, B] = Priv_BaseFcnParam(Xc,iXi);
        end
        
        function [GP,GW] = IntegrationScheme(order)
            % [GP,GW] = IntegrationScheme(order)
            % Gauss Points and Weights on the domain [0,1]
            
            [x,w] = gauss(order, 0,1);
            %x in [0,1]
            [X,Y,Z]=meshgrid(x,x,x);
            GP = [X(:),Y(:),Z(:)];
            
            [WX,WY,WZ]=meshgrid(w,w,w);
            GW = [WX(:),WY(:),WZ(:)];
            
            GW = prod(GW,2);
        end
        
    end
 
end

function [fi, detJ, B] = Priv_BaseFcnParam(XC,iXi)
% [fi, detJ, B] = BaseFcnParam(iel,iXi);
% Unit cube with node one in 0,0,0
%     8-----7
%    /|    /|
%   5-----6 |
%   | 4...|.3
%   |/    |/
%   1-----2
%
% x = XC(:,1);
% y = XC(:,2);
% z = XC(:,3);
% 
% xi = iXi(1,:);
% eta = iXi(2,:);
% zeta = iXi(3,:);
%
% xi in [0,1]
% eta in [0,1]
% zeta in [0,1]

xi = iXi(1);
eta = iXi(2);
zeta = iXi(3);

X0 = XC(:,1);
Y0 = XC(:,2);
Z0 = XC(:,3);

fi = [-(eta - 1)*(xi - 1)*(zeta - 1);...
    xi*(eta - 1)*(zeta - 1);...
    -eta*xi*(zeta - 1);...
    eta*(xi - 1)*(zeta - 1);...
    zeta*(eta - 1)*(xi - 1);...
    -xi*zeta*(eta - 1);...
    eta*xi*zeta;...
    -eta*zeta*(xi - 1)].';


dfidxi = [-(eta - 1)*(zeta - 1);
    (eta - 1)*(zeta - 1);
    -eta*(zeta - 1);
    eta*(zeta - 1);
    zeta*(eta - 1);
    -zeta*(eta - 1);
    eta*zeta;
    -eta*zeta];

dfideta = [-(xi - 1)*(zeta - 1);
    xi*(zeta - 1);
    -xi*(zeta - 1);
    (xi - 1)*(zeta - 1);
    zeta*(xi - 1);
    -xi*zeta;
    xi*zeta;
    -zeta*(xi - 1)];

dfidzeta = [-(eta - 1)*(xi - 1);
    xi*(eta - 1);
    -eta*xi;
    eta*(xi - 1);
    (eta - 1)*(xi - 1);
    -xi*(eta - 1);
    eta*xi;
    -eta*(xi - 1)];


J = [dfidxi.'*X0(:), dfidxi.'*Y0(:), dfidxi.'*Z0(:);...
    dfideta.'*X0(:), dfideta.'*Y0(:), dfideta.'*Z0(:);...
    dfidzeta.'*X0(:), dfidzeta.'*Y0(:), dfidzeta.'*Z0(:)];

detJ = det(J);

B = J\[dfidxi.';dfideta.';dfidzeta.'];

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
