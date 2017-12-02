% Testing code!
%
% Material_1 and Material_3 works as intended
% Material_2 needs more testing!

clear

dstrain = [0.001,0,0,0,0,0];

properties = [210000, 0.3];
[cauchy_stress_tensor,constitutive_tensor] = Material_1(dstrain,properties);

properties = [210000, 0.3, 1, 120, 7500];
[cauchy_stress_tensor,constitutive_tensor] = Material_2(dstrain,properties);

properties = [282692.307692308 121153.846153846 121153.846153846 0 0 0 282692.307692308 121153.846153846 0 0 0 282692.307692308 0 0 0 80769.2307692308 0 0 80769.2307692308 0 80769.2307692308];
[cauchy_stress_tensor,constitutive_tensor] = Material_3(dstrain,properties);