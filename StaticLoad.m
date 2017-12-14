function f = StaticLoad(model,step)
    % f = StaticLoad(model,step)

    mesh = model.mesh;
    dofs = size(mesh.P,2);
    neq = mesh.nnod*dofs;
    f = zeros(neq,1);
    nLoads = length(mesh.Loads);
    
    switch class(mesh)
        case 'C3D8_Mesh' 
            EdgeMesh.IntegrationScheme = @T3D2_Mesh.IntegrationScheme;
            ContinuumMesh.IntegrationScheme = @C3D8_Mesh.IntegrationScheme;
            SurfaceMesh.IntegrationScheme = @M3D4_Mesh.IntegrationScheme;
            
            EdgeMesh.BaseFcnParam_Static = @EdgeMesh.BaseFcnParam_Static2;
            ContinuumMesh.BaseFcnParam_Static = @C3D8_Mesh.BaseFcnParam_Static;
            SurfaceMesh.BaseFcnParam_Static = @M3D4_Mesh.BaseFcnParam_Static;
        otherwise
            error('Mesh type not implemented!')
    end
    
    
    
    if dofs==3
        %% 3D
        for iLoadType = 1:nLoads
            iLoad = mesh.Loads(iLoadType);
            nodes = iLoad.Set;
            [nele,knod] = size(nodes);
            switch iLoad.Type
                case LoadType.EdgeLoad
                    % User defined depending on Continuum mesh type.
                    [GP,GW] = EdgeMesh.IntegrationScheme(iLoad.IntegrationOrder);
                    BaseFunction = EdgeMesh.BaseFcnParam_Static;
                case LoadType.PointLoad
                    BaseFunction = @BaseFcnOneNode;
                case LoadType.BodyLoad
                    [GP,GW] = ContinuumMesh.IntegrationScheme(iLoad.IntegrationOrder);
                    BaseFunction = ContinuumMesh.BaseFcnParam_Static;
                case LoadType.SurfaceLoad
                    [GP,GW] = SurfaceMesh.IntegrationScheme(iLoad.IntegrationOrder);
                    BaseFunction = SurfaceMesh.BaseFcnParam_Static;
                otherwise
                    error('Load type not implemented')
            end
            
            
            %% Uses user defined base functions
            dir = iLoad.Direction;
            fval = iLoad.Magnitude;
            
            Bw = zeros(3,knod*dofs);
            
            for iel = 1:nele
                inods = nodes(iel,:);
                ieqs = zeros(3*length(inods),1);
                ieqs(1:3:end) = inods*3-2;
                ieqs(2:3:end) = inods*3-1;
                ieqs(3:3:end) = inods*3-0;
                
                Xc = mesh.P(inods,:);
                fe = zeros(knod*dofs,1);
                for ig = 1:size(GP,1)
                    iXi = GP(ig,:); iW = GW(ig);
                    [fi, detJ] = BaseFunction(Xc,iXi);
                    
                    Bw(1,1:3:end)=fi;
                    Bw(2,2:3:end)=fi;
                    Bw(3,3:3:end)=fi;
                    fe = fe + Bw'*dir(:)*fval*detJ*iW;
                end
                f(ieqs) = f(ieqs) + fe;
                
            end
            
        end
        
        
%         quiver3(mesh.P(:,1),mesh.P(:,2),mesh.P(:,3),F(1:3:end),F(2:3:end),F(3:3:end),0)
        
    elseif dofs==2
        %% 2D
        error('2D code is not implemented!')
    elseif dofs == 1
        %% 1D
        error('1D code is not implemented!')
    else
        error([num2str(dofs),'-dimensional problems not implemented!'])
    end
    
    

end

function [fi,detJ] = BaseFcnOneNode(Xc,iXi)
    fi = 1;
    detJ = 1;
end