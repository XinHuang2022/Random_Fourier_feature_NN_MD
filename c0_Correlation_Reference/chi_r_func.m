    
	function chi_value = chi_r_func( r, R_a, alpha_chi )
		
		%{
		if( r > R_b )
			chi_value = 0;
		elseif( r < R_a )
			chi_value = 1;
		else
			A_mat = [ 3 * R_b^2, 2 * R_b, 1, 0; 3 * R_a^2, 2 * R_a, 1, 0; R_b^3, R_b^2, R_b, 1; R_a^3, R_a^2, R_a, 1 ];
			b_vec = [ 0; 0; 0; 1 ];
			c_para = A_mat \ b_vec;
			a = c_para( 1, 1 );  b = c_para( 2, 1 );  c = c_para( 3, 1 );  d = c_para( 4, 1 );
			
			chi_value = a * r^3 + b * r^2 + c * r + d;
		end
		%}
		%{
		if( r < R_a )
			chi_value = 1;
		else
			chi_value = exp( -( r - R_a )^2 / ( 2 * alpha_chi ) );
		end
		%}
		chi_value = ( r < R_a ) + ( r >= R_a ) .* ( exp( -( r - R_a ).^2 / ( 2 * alpha_chi ) ) );
		
	end
	
	
	
	
	
	
	
	