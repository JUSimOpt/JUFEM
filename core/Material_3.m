% Anisotropic Elasticity, rate- and temperature independent for small deformations
%
% Inputs: 
%  stress [s11 s22 s33 s12 s13 s23]
%  dstrain [e11 e22 e33 e12 e13 e23]
%  properties [C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66] = [D1111 D1122 D1133 D1112 D1113 D1123 D2222 D2233 D2212 D2213 D2223 D3333 D3312 D3313 D3323 D1212 D1213 D1223 D1313 D1323 D2323]
%
% Outputs: 
%  stress [s11 s22 s33 s12 s13 s23]
%  constitutive_tensor [C11 C12 C13 C14 C15 C16] = [D1111 D1122 D1133 D1112 D1113 D1123] 
%                      [C21 C22 C23 C24 C25 C26]   [D2211 D2222 D2233 D2212 D2213 D2223]
%                      [C31 C32 C33 C34 C35 C36]   [D3311 D3322 D3333 D3312 D3313 D3323]
%                      [C41 C42 C43 C44 C45 C46]   [D1211 D1222 D1233 D1212 D1213 D1223]
%                      [C51 C52 C53 C54 C55 C56]   [D1311 D1322 D1333 D1312 D1313 D1323]
%                      [C61 C62 C63 C64 C65 C66]   [D2311 D2322 D2333 D2312 D2313 D2323]

function [stress,constitutive_tensor] = Material_3(stress,dstrain,properties)

    constitutive_tensor = [properties(1) properties(2) properties(3) properties(4) properties(5) properties(6);
                           properties(2) properties(7) properties(8) properties(9) properties(10) properties(11);
                           properties(3) properties(8) properties(12) properties(13) properties(14) properties(15);
                           properties(4) properties(9) properties(13) properties(16) properties(17) properties(18);
                           properties(5) properties(10) properties(14) properties(17) properties(19) properties(20);
                           properties(6) properties(11) properties(15) properties(18) properties(20) properties(21)];
    
    for i = 4:6
        dstrain(i) = dstrain(i)*2;
    end

	stress = stress + constitutive_tensor * dstrain;

end