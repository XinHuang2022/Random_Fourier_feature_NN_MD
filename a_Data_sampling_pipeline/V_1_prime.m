    
	
	function V_1_gradient = V_1_prime( x, alpha, gamma, R_a, R_c, sigma_a, sigma_c )
	    
		x_1 = x( :, 1 );
		x_2 = x( :, 2 );
		% x_3 = x( 1, 3 );
		r = sqrt( x_1.^2 + x_2.^2 );
		
		V_value = V_func( x, alpha, gamma, R_a, sigma_a );
		V_grad = V_prime( x, alpha, gamma, R_a, sigma_a );
		chi_1_value = chi_r_func( r, R_c, sigma_c );
		dchi_1_dr_value = dchi_dr_func( r, R_c, sigma_c );
	
		dV_dx1 = V_grad( :, 1 );
		dV_dx2 = V_grad( :, 2 );
		dV_1_dx1 = dV_dx1 .* chi_1_value + V_value .* dchi_1_dr_value .* x_1 ./ r;
		dV_1_dx2 = dV_dx2 .* chi_1_value + V_value .* dchi_1_dr_value .* x_2 ./ r;
		
		V_1_gradient = [ dV_1_dx1, dV_1_dx2 ];
		
	end
		
		
		
	
	
	
	
	
	
	
	