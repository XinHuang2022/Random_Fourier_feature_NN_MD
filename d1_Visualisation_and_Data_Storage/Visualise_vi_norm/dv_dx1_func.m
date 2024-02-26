    
	function dv_dx1_value = dv_dx1_func( x1, x2, alpha, gamma, R_a, sigma_a )
		
		R = sqrt( x1.^2 + x2.^2 );
		dv_dx1_value = x1 + gamma * ( x2 .* cos( x1 .* x2 ) .* chi_r_func( R, R_a, sigma_a ) + sin( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) .* x1 ./ R );
		
	end	
		
		
	
	
	
	
	
	
	
	
	