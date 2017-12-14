function [Ke,ge] = ElementData1PointIntegration_3D(materialData,Xc,ue,ieqs,GP,GW,baseFcnParam,iTime,nTimeIncrements,delTime)

    [fi, B, Bxi, Beta, Bzeta, detJ] = baseFcnParam(Xc);
    
    % B matrices and their xi-derivatives
    [BN,BL,BNXi,BNEta,BNZeta,BLXi,BLEta,BLZeta] = ComputeB_1Point(B,Bxi,Beta,Bzeta,ue);
    
    gradU = BN*ue;
    
    [D_iso,Sv_iso,D_vol,Sv_vol]=materialData.materialFcn(gradU);
    
    D = D_iso + D_vol; % Constitutive tensor on Voigt form
    
    % Stress tensor
    ind = [1,4,6;
           4,2,5;
           6,5,3];
    S_iso = Sv_iso(ind);
    S_vol = Sv_vol(ind);
    Sv = Sv_iso + Sv_vol;
    
    T_iso = [S_iso,zeros(3),zeros(3);
             zeros(3),S_iso,zeros(3);
             zeros(3),zeros(3),S_iso];
    T_vol = [S_vol,zeros(3),zeros(3);
             zeros(3),S_vol,zeros(3);
             zeros(3),zeros(3),S_vol];
         
    T = T_iso + T_vol;
    
    Stab_iso = 1/12*(BLXi'*D_iso*BLXi + BNXi'*T_iso*BNXi + ...
                   BLEta'*D_iso*BLEta + BNEta'*T_iso*BNEta + ...
                   BLZeta'*D_iso*BLZeta + BNZeta'*T_iso*BNZeta);          


    Ke =  detJ*( (BL'*D*BL+BN'*T*BN) + Stab_iso);
    ge =  detJ*( (BL'*Sv) + Stab_iso*ue);
               
    

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