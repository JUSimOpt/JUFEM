function [cauchy_stress_tensor,constitutive_tensor] = Material_1(dstrain,properties)

	lam = properties(1) * properties(2) / ((1 - 2 * properties(2)) * (1 + properties(2)))
	mu = properties(1) / (2 * (1 + properties(2)))
	
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

end