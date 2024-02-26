	function dv_dx2_value = dv_dx2_func( x1, x2, alpha, gamma, R_a, sigma_a )
		
		R = sqrt( x1.^2 + x2.^2 );
		dv_dx2_value = alpha * x2 + gamma * ( x1 .* cos( x1 .* x2 ) .* chi_r_func( R, R_a, sigma_a ) + sin( x1 .* x2 ) .* dchi_dr_func( R, R_a, sigma_a ) .* x2 ./ R );
		
	end	