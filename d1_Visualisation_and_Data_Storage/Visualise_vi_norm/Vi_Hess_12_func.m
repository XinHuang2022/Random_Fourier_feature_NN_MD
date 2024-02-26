    
	function Vi_Hess_12_value = Vi_Hess_12_func( x1, x2, alpha, gamma, R_a, sigma_a, R_b, sigma_b )
		
		R = sqrt( x1.^2 + x2.^2 );
		Vi_Hess_12_value = dchi_dr_func( R, R_b, sigma_b ) .* ( x1 ./ R .* dv_dx2_func( x1, x2, alpha, gamma, R_a, sigma_a ) + x2 ./ R .* dv_dx1_func( x1, x2, alpha, gamma, R_a, sigma_a ) )...
						   - V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a ) .* dchi_dr_func( R, R_b, sigma_b ) .* ( x1 .* x2 ./ R.^3 )...
						   + chi_r_func( R, R_b, sigma_b ) .* d2v_dx1_dx2_func( x1, x2, alpha, gamma, R_a, sigma_a )...
						   + V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a ) .* d2chi_dr_2_func( R, R_b, sigma_b ) .* x1 .* x2 ./ R.^2;
		
		
		
		
	
	
	
	
	