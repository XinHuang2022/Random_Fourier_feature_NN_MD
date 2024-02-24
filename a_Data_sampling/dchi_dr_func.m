    
	function dchi_dr = dchi_dr_func( r, R, alpha_chi )
		
		%{
		if( ( r > R_b ) || ( r < R_a ) )
			dchi_dr = 0;
		else
			A_mat = [ 3 * R_b^2, 2 * R_b, 1, 0; 3 * R_a^2, 2 * R_a, 1, 0; R_b^3, R_b^2, R_b, 1; R_a^3, R_a^2, R_a, 1 ];
			b_vec = [ 0; 0; 0; 1 ];
			c_para = A_mat \ b_vec;
			a = c_para( 1, 1 );  b = c_para( 2, 1 );  c = c_para( 3, 1 );  d = c_para( 4, 1 );
			
			dchi_dr = 3 * a * r^2 + 2 * b * r + c;
		end
		%}
		%{
		if( r < R_a )
			dchi_dr = 0;
		else
			dchi_dr = -( r - R_a ) / alpha_chi * exp( -( r - R_a )^2 / ( 2 * alpha_chi ) );
		end
		%}
		dchi_dr = ( r > R ) .* ( -( r - R ) / alpha_chi .* exp( -( r - R ).^2 / ( 2 * alpha_chi ) ) );
		
    end
		
		
	
	
	
	
	
	
	
	
	
	