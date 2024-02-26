    function d2chi_dr_2_value = d2chi_dr_2_func( R, R_ref, sigma_ref )
		
		d2chi_dr_2_value = exp( -( R - R_ref ).^2 / ( 2 * sigma_ref ) ) .* ( ( R - R_ref ).^2 / ( sigma_ref^2 ) - 1 / sigma_ref );
		
	end
		
		
		
		
		
		
		
		
		