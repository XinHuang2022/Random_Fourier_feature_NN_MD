    
	function d2v_dx1_dx2_value = d2v_dx1_dx2_func( x1, x2, alpha, gamma, R_a, sigma_a )
	
		R = sqrt( x1.^2 + x2.^2 );
		d2v_dx1_dx2_value = gamma * ( ( cos( x1 .* x2 ) + x1 .* x2 .* sin( x1 .* x2 ) ) .* chi_r_func( R, R_a, sigma_a )...
									  + cos( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) .* R...
									  + sin( x1 .* x2 ) .* d2chi_dr_2_func( R, R_a, sigma_a ) .* x1 .* x2 ./ R.^2 ...
									  - sin( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) .* x1 .* x2 ./ R.^3 );
									  
	end
									
									
									
		
		
		
		
		
	
	
	