function [cauchy_stress_tensor,constitutive_tensor,history_variables] = Material_2(dstrain,properties,history_variables)

	lam = properties(1) * properties(2) / ((1 - 2 * properties(2)) * (1 + properties(2)))
	mu = properties(1) / (2 * (1 + properties(2)))
	mixed = 1 % Mixed hardening?
	
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
			cauchy_stress_tensor(i) = cauchy_stress_tensor(i) + constitutive_tensor(i,j) * dstran(j);
		end
	end
     
% back stress
	for i = 1:6
		alpha(i) = history_variables(i)
	end

% trial stress
	for i = 1:6
		trial_stress(i) = cauchy_stress_tensor(i) - alpha(i)
	end
	hydrostatic_stress = (1/3) * (cauchy_stress_tensor(1) + cauchy_stress_tensor(2) + cauchy_stress_tensor(3))
	for i = 1:3
		trial_stress(i) = trial_stress(i) - hydrostatic_stress
	end
	effective_trial_stress = sqrt(trial_stress(1)^2 + trial_stress(2)^2 + trial_stress(3)^2 + 2*trial_stress(4)^2 + 2*trial_stress(5)^2 + 2*trial_stress(6)^2)
	for i = 1:6
		unit_trial_stress(i) = trial_stress(i) / effective_trial_stress
	end

% plastic corrector
	gammadt_n = history_variables(7)
	if (abs(effective_trial_stress) > sqrt(2/3) * fkappa(properties,history_variables(7))) then
		while abs(g) > tol
			junk = history_variables(7) + sqrt(2/3) * gammadt
			g = effective_trial_stress - sqrt(2/3) * fkappa(properties,junk) - sqrt(2/3) * (H(properties,junk) - H(properties,gammadt_n)) - 2 * mu * gammadt
			dg = -2 * mu * (1 + ((dH(properties,junk) + dkappa(properties,junk)) / (3*mu)))
			gammadt = gammadt - (g / dg)
		end
		history_variables(7) = history_variables(7) + sqrt(2/3) * gammadt

		for i = 1:6
			alpha(i) = alpha(i) + sqrt(2/3) * (H(properties,history_variables(7))-H(properties,gammadt_n))*unit_trial_stress(i)
			cauchy_stress_tensor(i) = sqrt(2/3) * fkappa(properties,history_variables(7)) * unit_trial_stress(i) + alpha(i)
			history_variables(i) = alpha(i)
		end
		for i = 1:3
			cauchy_stress_tensor(i) = cauchy_stress_tensor(i) + hydrostatic_stress
		end

% consistent tangent
		for i = 1:6
			for j = 1:6
				vi4(i,j) = 0 % Prevents NaN entries in the stiffness matrix (fortran ftw!?)
			end
		end
		for i = 1:3
			vi2(i) = 1
			vi4(i,i) = 1
		end
		for i = 4:6
			vi2(i) = 0
			vi4(i,i) = 1/2
		end
		for i = 1:6
			for j = 1:6
				vi2prod(i,j) = vi2(i) * vi2(j)
			end
		end
		for i = 1:6
			for j = 1:6
				trprod(i,j) = unit_trial_stress(i) * unit_trial_stress(j)
			end
		end
		beta = sqrt(2/3) * ((fkappa(properties,history_variables(7)) + (H(properties,history_variables(7)) - H(properties,gammadt_n))) / effective_trial_stress)
		gamma = (1 / (1 + ((dkappa(properties,history_variables(7)) + dH(properties,history_variables(7))) / (3 * mu)))) - (1 - beta)
		for i = 1:6
			for j = 1:6
				ddsdde(i,j) = (lam + (2/3) * mu) * vi2prod(i,j) + 2 * mu * beta * (vi4(i,j) - (1/3) * vi2prod(i,j)) - 2 * mu * gamma * trprod(i,j)
			end
		end
	end

% HARDENING RULE
	function fkappa(properties,peeq)
		fkappa = properties(3) + mixed * properties(4) * peeq
	end
	function dkappa(properties,peeq)
		dkappa = mixed * properties(4) * peeq
	end

% PLASTIC MODULUS
	function H(properties,peeq)
		H = (1 - mixed) * properties(4) * peeq
	end
	function dH(properties,peeq)
		dH = (1 - mixed) * 0
	end
      
end
