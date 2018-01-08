function [Ke,ge,histvar] = UserElement_3D_Full(histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter,baseFcnParam)

%     loadFrac = iTime/nTimeIncrements;

    nLocDofs = length(ieqs);
    dofs = size(GP,2);
    ge = zeros(nLocDofs,1); % element internal load vector
%     re = fe; % element resultant load vector
    Ke = zeros(nLocDofs,nLocDofs); % element tangential stiffness matrix
%     Bw = zeros(3,nLocDofs);
    
    % External load
%     fe = fe*loadFrac;
    ind = [1,4,6;
           4,2,5;
           6,5,3];
    I = eye(3);
    Iv = I(:);
    
    E = zeros(3); Stress = E;
    
    nIntegrationPoints = size(GP,1);
    for iG = 1:nIntegrationPoints
        iXi = GP(iG,:); iW = GW(iG);
        [fi, detJ, B] = baseFcnParam(Xc,iXi);
%         [BN,BL,BNXi,BNEta,BNZeta,BLXi,BLEta,BLZeta] = ComputeB1Point(B,Bxi,Beta,Bzeta,Ue);
        [BN,BL] = ComputeB(B,ue);
        gradU = BN*ue;
        
        
        F=reshape(Iv+gradU,3,3)';
        C=F'*F;
        E = (C-I)/2; %Engineering strain
        
        
%         [D_iso,D_vol,Sv_iso,Sv_vol]=NeoHook3D_2PK(H,C10,D1)
        [D_iso,Sv_iso,D_vol,Sv_vol]=materialData.materialFcn(gradU);
        
        D = D_iso + D_vol; % Constitutive tensor on Voigt form
        Sv = Sv_iso + Sv_vol; % Stress tensor on Voigt form
        
        % Stress tensor 
        S = Sv(ind);
        T = [S,zeros(3),zeros(3);
             zeros(3),S,zeros(3);
             zeros(3),zeros(3),S];
        
        % Integrate the tangential matrix 
        Ke = Ke + (BL'*D*BL + BN'*T*BN) * detJ*iW;
        
        % Integrate the internal forces
        ge = ge + (BL'*Sv) * detJ*iW;
        
    end
%     re = (fe-ge);
1;
end

function [BN,BL] = ComputeB(B,Ue)

    fix = B(1,:);
    fiy = B(2,:);
    fiz = B(3,:);
    BN = zeros(9,3*8);
    BN(1,1:3:end)=fix;
    BN(2,1:3:end)=fiy;
    BN(3,1:3:end)=fiz;
    BN(4,2:3:end)=fix;
    BN(5,2:3:end)=fiy;
    BN(6,2:3:end)=fiz;
    BN(7,3:3:end)=fix;
    BN(8,3:3:end)=fiy;
    BN(9,3:3:end)=fiz;
    H=BN*Ue;
    del=[1 0 0 0 1 0 0 0 1]';
    F=reshape(del+H,3,3)';

    
    
    fix = B(1,:);
    fiy = B(2,:);
    fiz = B(3,:);
    BL=zeros(6,3*8);
    BL(1,1:3:end)=F(1,1)*fix;
    BL(1,2:3:end)=F(2,1)*fix;
    BL(1,3:3:end)=F(3,1)*fix;
    BL(2,1:3:end)=F(1,2)*fiy;
    BL(2,2:3:end)=F(2,2)*fiy;
    BL(2,3:3:end)=F(3,2)*fiy;
    BL(3,1:3:end)=F(1,3)*fiz;
    BL(3,2:3:end)=F(2,3)*fiz;
    BL(3,3:3:end)=F(3,3)*fiz;
    BL(4,1:3:end)=F(1,1)*fiy+F(1,2)*fix;
    BL(4,2:3:end)=F(2,1)*fiy+F(2,2)*fix;
    BL(4,3:3:end)=F(3,1)*fiy+F(3,2)*fix;
    BL(5,1:3:end)=F(1,2)*fiz+F(1,3)*fiy;
    BL(5,2:3:end)=F(2,2)*fiz+F(2,3)*fiy;
    BL(5,3:3:end)=F(3,2)*fiz+F(3,3)*fiy;
    BL(6,1:3:end)=F(1,3)*fix+F(1,1)*fiz;
    BL(6,2:3:end)=F(2,3)*fix+F(2,1)*fiz;
    BL(6,3:3:end)=F(3,3)*fix+F(3,1)*fiz;
    
    
end

