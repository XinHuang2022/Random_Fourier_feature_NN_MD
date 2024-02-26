    
	function chi_r_value = chi_r_func( R, R_a, sigma_a )
		
		chi_r_value = ( R <= R_a ) + ( R > R_a ) .* exp( -( R - R_a ).^2 / ( 2 * sigma_a ) );
		
	end
		
		
		
		
	
	
	
	
	