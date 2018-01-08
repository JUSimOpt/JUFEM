function [Ke,ge,histvar] = UserElement_3D_Reduced(histvar,materialData,Xc,ue,ieqs,GP,GW,cInc,time,dTime,iter,baseFcnParam)


    I = eye(3);
    Iv = I(:);  
    
    nLocDofs = length(ieqs);
    dofs = size(GP,2);
    ge = zeros(nLocDofs,1); % element internal load vector
    Ke = zeros(nLocDofs,nLocDofs); % element tangential stiffness matrix
    ge_iso = ge;
    Ke_iso = Ke;
    Ke_vol = Ke_iso;
    ge_vol = ge_iso;
    ind = [1,4,6;
           4,2,5;
           6,5,3];
    nGP = size(GP,1);
    
      
    for i = 1:nGP
       iXi =  GP(i,:); iw = GW(i);
       [fi, detJ, B] = baseFcnParam(Xc,iXi);
       [BN,BL] = ComputeB(B,ue);
        gradU = BN*ue;
        [D_iso,Sv_iso,D_vol,Sv_vol]=materialData.materialFcn(gradU);
        % Stress tensor
        
        S_iso = Sv_iso(ind);
        S_vol = Sv_vol(ind);
%         
        T_iso = [S_iso,zeros(3),zeros(3);
            zeros(3),S_iso,zeros(3);
            zeros(3),zeros(3),S_iso];
        T_vol = [S_vol,zeros(3),zeros(3);
            zeros(3),S_vol,zeros(3);
            zeros(3),zeros(3),S_vol];
%         
        D = D_iso + D_vol; % Constitutive tensor on Voigt form
        Sv = Sv_iso + Sv_vol; % Stress tensor on Voigt form
        
        % Stress tensor 
        S = Sv(ind);
        
        F=reshape(Iv+gradU,3,3)'; %Deformation gradient
        C=F'*F; %Right Cauchy-Green tensor
        E = (C-I)/2; %Engineering strain
        
        histvar.ip(i).S = S;
        histvar.ip(i).EE = E;
        
%         T = [S,zeros(3),zeros(3);
%              zeros(3),S,zeros(3);
%              zeros(3),zeros(3),S];
        
        Ke_iso = Ke_iso + (BL'*D_iso*BL + BN'*T_iso*BN) * detJ*iw;
        Ke_vol = Ke_vol + (BL'*D_vol*BL + BN'*T_vol*BN) * detJ*iw;
        
        ge_iso = ge_iso + (BL'*Sv_iso) * detJ*iw;
        ge_vol = ge_vol + (BL'*Sv_vol) * detJ*iw;
        
        
    end
    1;
    % Midpoint
    iXi = [1,1,1]/2; iw = 1;
    [fi, detJ, B] = baseFcnParam(Xc,iXi);
    [BN,BL] = ComputeB(B,ue);
    gradU = BN*ue;
    [D_iso,Sv_iso,D_vol,Sv_vol]=materialData.materialFcn(gradU);
    S_vol = Sv_vol(ind);
    T_vol = [S_vol,zeros(3),zeros(3);
        zeros(3),S_vol,zeros(3);
        zeros(3),zeros(3),S_vol];
    
    Ke_vol = (BL'*D_vol*BL + BN'*T_vol*BN) * detJ*iw;
    ge_vol = (BL'*Sv_vol) * detJ*iw;
    
    
    Ke = Ke_iso + Ke_vol;
    ge = ge_iso + ge_vol;
    
    
    
    
%     [fi, B, Bxi, Beta, Bzeta, detJ] = baseFcnParam_1P(Xc);
%     
% %     [fi, detJ, B] = baseFcnParam(Xc,iXi);
%     
%     % B matrices and their xi-derivatives
%     [BN,BL,BNXi,BNEta,BNZeta,BLXi,BLEta,BLZeta] = ComputeB_1Point(B,Bxi,Beta,Bzeta,ue);
%     
%     gradU = BN*ue;
%     
%     [D_iso,Sv_iso,D_vol,Sv_vol]=materialData.materialFcn(gradU);
%     
%     D = D_iso + D_vol; % Constitutive tensor on Voigt form
%     
%     % Stress tensor
%     ind = [1,4,6;
%            4,2,5;
%            6,5,3];
%     S_iso = Sv_iso(ind);
%     S_vol = Sv_vol(ind);
%     Sv = Sv_iso + Sv_vol;
%     
%     T_iso = [S_iso,zeros(3),zeros(3);
%              zeros(3),S_iso,zeros(3);
%              zeros(3),zeros(3),S_iso];
%     T_vol = [S_vol,zeros(3),zeros(3);
%              zeros(3),S_vol,zeros(3);
%              zeros(3),zeros(3),S_vol];
%          
%     T = T_iso + T_vol;
%     
%     Stab_iso = 1/12*(BLXi'*D_iso*BLXi + BNXi'*T_iso*BNXi + ...
%                    BLEta'*D_iso*BLEta + BNEta'*T_iso*BNEta + ...
%                    BLZeta'*D_iso*BLZeta + BNZeta'*T_iso*BNZeta);    
%                
%     Stab_vol = 1/12*(BLXi'*D_vol*BLXi + BNXi'*T_vol*BNXi + ...
%                    BLEta'*D_vol*BLEta + BNEta'*T_vol*BNEta + ...
%                    BLZeta'*D_vol*BLZeta + BNZeta'*T_vol*BNZeta);   
% 
%     Ke =  detJ*( (BL'*D*BL+BN'*T*BN) + Stab_iso );
%     ge =  detJ*( (BL'*Sv) + (Stab_iso)*ue);
               
    

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



function [BN,BL,BNXi,BNEta,BNZeta,BLXi,BLEta,BLZeta] = ComputeB_1Point(B,Bxi,Beta,Bzeta,Ue)

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

    fix = Bxi(1,:);
    fiy = Bxi(2,:);
    fiz = Bxi(3,:);
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
    BLXi = BL;

    fix = Beta(1,:);
    fiy = Beta(2,:);
    fiz = Beta(3,:);
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
    BLEta = BL;

    fix = Bzeta(1,:);
    fiy = Bzeta(2,:);
    fiz = Bzeta(3,:);
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
    BLZeta = BL;
    
    
    fix = Bxi(1,:);
    fiy = Bxi(2,:);
    fiz = Bxi(3,:);
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
    BNXi = BN;
    
    fix = Beta(1,:);
    fiy = Beta(2,:);
    fiz = Beta(3,:);
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
    BNEta = BN;
    
    fix = Bzeta(1,:);
    fiy = Bzeta(2,:);
    fiz = Bzeta(3,:);
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
    BNZeta = BN;
    
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
    
end