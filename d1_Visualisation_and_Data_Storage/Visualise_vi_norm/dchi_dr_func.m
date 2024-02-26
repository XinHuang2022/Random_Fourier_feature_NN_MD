	
	function dchi_dr_value = dchi_dr_func( R, R_cf, sigma_cf )
		
		dchi_dr_value = ( R > R_cf ) .* ( -( R - R_cf ) .* exp( -( R - R_cf ).^2 / ( 2 * sigma_cf ) ) / sigma_cf );
	
	end
		
		
		
		
	
	
	
	
	
	
	
	
	
	