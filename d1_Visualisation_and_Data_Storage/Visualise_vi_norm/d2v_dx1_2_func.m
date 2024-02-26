    
	function d2v_dx1_2_value = d2v_dx1_2_func( x1, x2, alpha, gamma, R_a, sigma_a )
		
		R = sqrt( x1.^2 + x2.^2 );
		d2v_dx1_2_value = 1 + gamma * ( x2 .* cos( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) .* x1 ./ R...
										- x2.^2 .* sin( x1 .* x2 ) .* chi_r_func( R, R_a, sigma_a )...
										+ ( x2 .* cos( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) + sin( x1 .* x2 ) .* d2chi_dr_2_func( R, R_a, sigma_a ) .* x1 ./ R ) .* ( x1 ./ R )...
										+ sin( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) .* x2.^2 ./ R.^3 );
									   
	
	end
	
	
	
	
	
	
	
	
	