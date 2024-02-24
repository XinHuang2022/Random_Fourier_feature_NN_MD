    
    function V_value = V_func( x, alpha, gamma, R_minus_2, alpha_chi )
        
        % V_value = x * x' / 2 ;
		x_1 = x( :, 1 );
		x_2 = x( :, 2 );
		% x_3 = x( 1, 3 );
		r = sqrt( x_1.^2 + x_2.^2 );
		
		% V_value = x_1^2 / 2 + alpha * x_2^2 / 2 + gamma * sin( x_1 * x_2 ) + x_3^2 / 2;
		V_value = x_1.^2 / 2 + alpha * x_2.^2 / 2 + gamma * sin( x_1 .* x_2 ) .* chi_r_func( r, R_minus_2, alpha_chi );
		
        
    end
    
    
    
    
    
    
    
    
    