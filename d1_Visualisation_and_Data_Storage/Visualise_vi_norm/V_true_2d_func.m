    
	function V_2d_input = V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a )
		
		% alpha = sqrt( 2 );
		% gamma = 2;
		
		%{
		A_mat = [ 3 * R_b^2, 2 * R_b, 1, 0; 3 * R_a^2, 2 * R_a, 1, 0; R_b^3, R_b^2, R_b, 1; R_a^3, R_a^2, R_a, 1 ];
		b_vec = [ 0; 0; 0; 1 ];
		c_para = A_mat \ b_vec;
		a = c_para( 1, 1 );  b = c_para( 2, 1 );  c = c_para( 3, 1 );  d = c_para( 4, 1 );
		%}
		
		R = sqrt( x1.^2 + x2.^2 );
		
		V_2d_input = ( x1.^2 / 2 + alpha * x2.^2 / 2 + gamma * sin( x1 .* x2 ) .* chi_r_func( R, R_a, sigma_a ) ) ;
		
		% V_2d_input = ( R < R_a ) .* ( x1.^2 / 2 + alpha * x2.^2 / 2 + gamma * sin( x1 .* x2 ) ) ...
					 % + ( R >= R_a ) .* ( x1.^2 / 2 + alpha * x2.^2 / 2 + gamma * sin( x1 .* x2 ) .* exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) );
					 % + ( ( R >= R_a ) .* ( R <= R_b ) ) .* ( x1.^2 / 2 + alpha * x2.^2 / 2 + gamma * sin( x1 .* x2 ) .* ( a * R.^3 + b * R.^2 + c * R + d ) ) ...
					 % + ( R > R_b ) .* ( x1.^2 / 2 + alpha * x2.^2 / 2 );
		
    end
		
		
		
	
	
	
	
	
	