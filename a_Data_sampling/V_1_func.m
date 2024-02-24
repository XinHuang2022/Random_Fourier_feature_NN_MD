    
	
	function V_1_value = V_1_func( x, alpha, gamma, R_a, R_c, sigma_a, sigma_c )
		
		x_1 = x( :, 1 );
		x_2 = x( :, 2 );
		% x_3 = x( 1, 3 );
		r = sqrt( x_1.^2 + x_2.^2 );
		
		V_value = V_func( x, alpha, gamma, R_a, sigma_a );
		chi_1_value = chi_r_func( r, R_c, sigma_c );
		V_1_value =  V_value .* chi_1_value;
		
	end
		
		
		
		
	
	
	
	
	
	
	
	
	