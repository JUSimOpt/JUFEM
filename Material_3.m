function [cauchy_stress_tensor,constitutive_tensor] = Material_3(dstrain,properties)

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
			cauchy_stress_tensor(i) = cauchy_stress_tensor(i) + constitutive_tensor(i,j) * dstran(j);
		end
	end

% test (remove this!)
    for i = 1:6
		cauchy_stress_tensor_test(i) = cauchy_stress_tensor_test(i) + constitutive_tensor(i,1) * dstran(1) + constitutive_tensor(i,2) * dstran(2) + constitutive_tensor(i,3) * dstran(3) + constitutive_tensor(i,4) * dstran(4) + constitutive_tensor(i,5) * dstran(5) + constitutive_tensor(i,6) * dstran(6)
	end
      
end