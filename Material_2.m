% Isotropic Mixed Hardening Plasticity, rate- and temperature independent for small deformations
%
% Inputs: 
%  dstrain [e11,e22,e33,e12,e13,e23]
%  properties [E,nu,p,Y,H]
%
% Outputs: 
%  cauchy_stress_tensor [s11,s22,s33,s12,s13,s23]
%  constitutive_tensor [D1111,D1122,D1133,D1112,D1113,D1123]
%                      [D2211,D2222,D2233,D2212,D2213,D2223]
%                      [D3311,D3322,D3333,D3312,D3313,D3323]
%                      [D1211,D1222,D1233,D1212,D1213,D1223]
%                      [D1311,D1322,D1333,D1312,D1313,D1323]
%                      [D2311,D2322,D2322,D2312,D2313,D2323]

function [cauchy_stress_tensor,constitutive_tensor] = Material_2(dstrain,properties)

	lam = properties(1) * properties(2) / ((1 - 2 * properties(2)) * (1 + properties(2)));
	mu = properties(1) / (2 * (1 + properties(2)));
    
	cauchy_stress_tensor = zeros(1,6);
	constitutive_tensor = zeros(6,6);
    alpha = zeros(1,6);
    trial_stress = zeros(1,6);
    unit_trial_stress = zeros(1,6);
    vi2 = zeros(1,6);
    vi4 = zeros(6,6);
    vi2prod = zeros(6,6);
    trprod = zeros(6,6);
    g = 0;
    gammadt = 0;
    tol = 0.0001;
  
    history_variables = zeros(1,7);
    
% elastic predictor
	for i = 1:3
		for j = 1:3
			constitutive_tensor(i,j) = lam;
		end
		constitutive_tensor(i,i) = 2 * mu + lam;
	end
	for i = 4:6
		constitutive_tensor(i,i) = mu;
	end
	for i = 1:6
		for j = 1:6
			cauchy_stress_tensor(i) = cauchy_stress_tensor(i) + constitutive_tensor(i,j) * dstrain(j);
		end
	end
     
% back stress
	for i = 1:6
		alpha(i) = history_variables(i);
	end

% trial stress
	for i = 1:6
		trial_stress(i) = cauchy_stress_tensor(i) - alpha(i);
	end
	hydrostatic_stress = (1/3) * (cauchy_stress_tensor(1) + cauchy_stress_tensor(2) + cauchy_stress_tensor(3));
	for i = 1:3
		trial_stress(i) = trial_stress(i) - hydrostatic_stress;
	end
	effective_trial_stress = sqrt(trial_stress(1)^2 + trial_stress(2)^2 + trial_stress(3)^2 + 2*trial_stress(4)^2 + 2*trial_stress(5)^2 + 2*trial_stress(6)^2);
	for i = 1:6
		unit_trial_stress(i) = trial_stress(i) / effective_trial_stress;
	end

% plastic corrector
	gammadt_n = history_variables(7);
	if abs(effective_trial_stress) > sqrt(2/3) * fkappa(properties(3), properties(4), properties(5), history_variables(7))
        junk = history_variables(7) + sqrt(2/3) * gammadt;
		g = effective_trial_stress - sqrt(2/3) * fkappa(properties(3),properties(4),properties(5),junk) - sqrt(2/3) * (H(properties(3),properties(5),junk) - H(properties(3),properties(5),gammadt_n)) - 2 * mu * gammadt;
		while abs(g) > tol
			g = effective_trial_stress - sqrt(2/3) * fkappa(properties(3),properties(4),properties(5),junk) - sqrt(2/3) * (H(properties(3),properties(5),junk) - H(properties(3),properties(5),gammadt_n)) - 2 * mu * gammadt;
			dg = -2 * mu * (1 + ((dH(properties(3),junk) + dkappa(properties(3),properties(5),junk)) / (3*mu)));
			gammadt = gammadt - (g / dg);
		end
		history_variables(7) = history_variables(7) + sqrt(2/3) * gammadt;

		for i = 1:6
			alpha(i) = alpha(i) + sqrt(2/3) * (H(properties(3),properties(5),history_variables(7))-H(properties(3),properties(5),gammadt_n))*unit_trial_stress(i);
			cauchy_stress_tensor(i) = sqrt(2/3) * fkappa(properties(3),properties(4),properties(5),history_variables(7)) * unit_trial_stress(i) + alpha(i);
			history_variables(i) = alpha(i);
		end
		for i = 1:3
			cauchy_stress_tensor(i) = cauchy_stress_tensor(i) + hydrostatic_stress;
		end

% consistent tangent
		for i = 1:3
			vi2(i) = 1;
			vi4(i,i) = 1;
		end
		for i = 4:6
			vi2(i) = 0;
			vi4(i,i) = 1/2;
		end
		for i = 1:6
			for j = 1:6
				vi2prod(i,j) = vi2(i) * vi2(j);
			end
		end
		for i = 1:6
			for j = 1:6
				trprod(i,j) = unit_trial_stress(i) * unit_trial_stress(j);
			end
		end
		beta = sqrt(2/3) * ((fkappa(properties(3),properties(4),properties(5),history_variables(7)) + (H(properties(3),properties(5),history_variables(7)) - H(properties(3),properties(5),gammadt_n))) / effective_trial_stress);
		gamma = (1 / (1 + ((dkappa(properties(3),properties(5),history_variables(7)) + dH(properties(3),history_variables(7))) / (3 * mu)))) - (1 - beta);
		for i = 1:6
			for j = 1:6
				constitutive_tensor(i,j) = (lam + (2/3) * mu) * vi2prod(i,j) + 2 * mu * beta * (vi4(i,j) - (1/3) * vi2prod(i,j)) - 2 * mu * gamma * trprod(i,j);
			end
        end
	end
      
end

% HARDENING RULE
function a = fkappa(p,Y,H,peeq)
	a = Y + p * H * peeq;
end
function b = dkappa(p,H,peeq)
	b = p * H * peeq;
end

% PLASTIC MODULUS
function c = H(p,H,peeq)
	c = (1 - p) * H * peeq;
end
function d = dH(p,peeq)
	d = (1 - p) * 0;
end
