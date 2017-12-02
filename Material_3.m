% Anisotropic Elasticity, rate- and temperature independent for small deformations
%
% Inputs: 
%  dstrain [e11,e22,e33,2*e12,2*e13,2*e23]
%  properties [C11,C12,C13,C14,C15,C16,C22,C23,C24,C25,C26,C33,C34,C35,C36,C44,C45,C46,C55,C56,C66] = [D1111,D1122,D1133,D1112,D1113,D1123,D2222,D2233,D2212,D2213,D2223,D3333,D3312,D3313,D3323,D1212,D1213,D1223,D1313,D1323,D2323]
%
% Outputs: 
%  cauchy_stress_tensor [s11,s22,s33,s12,s13,s23]
%  constitutive_tensor [C11,C12,C13,C14,C15,C16] = [D1111,D1122,D1133,D1112,D1113,D1123] 
%                      [C21,C22,C23,C24,C25,C26]   [D2211,D2222,D2233,D2212,D2213,D2223]
%                      [C31,C32,C33,C34,C35,C36]   [D3311,D3322,D3333,D3312,D3313,D3323]
%                      [C41,C42,C43,C44,C45,C46]   [D1211,D1222,D1233,D1212,D1213,D1223]
%                      [C51,C52,C53,C54,C55,C56]   [D1311,D1322,D1333,D1312,D1313,D1323]
%                      [C61,C62,C63,C64,C65,C66]   [D2311,D2322,D2333,D2312,D2313,D2323]

function [cauchy_stress_tensor,constitutive_tensor] = Material_3(dstrain,properties)

	cauchy_stress_tensor = zeros(1,6);
	constitutive_tensor = zeros(6,6);

	for i = 1:6
		constitutive_tensor(1,i) = properties(i);
        constitutive_tensor(i,1) = properties(i);
	end
	for i = 1:5
		constitutive_tensor(2,i+1) = properties(i+6);
		constitutive_tensor(i+1,2) = properties(i+6);
	end
	for i = 1:4
		constitutive_tensor(3,i+2) = properties(i+11);
		constitutive_tensor(i+2,3) = properties(i+11);
	end
	for i = 1:3
		constitutive_tensor(4,i+3) = properties(i+15);
		constitutive_tensor(i+3,4) = properties(i+15);
	end
	for i = 1:2
		constitutive_tensor(5,i+4) = properties(i+18);
		constitutive_tensor(i+4,5) = properties(i+18);
	end
	constitutive_tensor(6,6) = properties(21);


	for i = 1:6
		for j = 1:6
			cauchy_stress_tensor(i) = cauchy_stress_tensor(i) + constitutive_tensor(i,j) * dstrain(j);
		end
	end

end