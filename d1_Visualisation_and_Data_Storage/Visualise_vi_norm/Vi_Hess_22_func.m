    
	    
	function Vi_Hess_22_value = Vi_Hess_22_func( x1, x2, alpha, gamma, R_a, sigma_a, R_b, sigma_b )
		
		R = sqrt( x1.^2 + x2.^2 );
		Vi_Hess_22_value = 2 * dchi_dr_func( R, R_b, sigma_b ) .* x2 ./ R .* dv_dx2_func( x1, x2, alpha, gamma, R_a, sigma_a )...
						   + chi_r_func( R, R_b, sigma_b ) .* d2v_dx2_2_func( x1, x2, alpha, gamma, R_a, sigma_a )...
						   + V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a ) .* d2chi_dr_2_func( R, R_b, sigma_b ) .* x2.^2 ./ R.^2 ...
						   + V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a ) .* dchi_dr_func( R, R_b, sigma_b ) .* x1.^2 ./ R.^3;
	
	end
		
		
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	