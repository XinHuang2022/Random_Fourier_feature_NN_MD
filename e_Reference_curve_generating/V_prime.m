    
    function V_gradient = V_prime( x, alpha, gamma, R_minus_2, alpha_chi )
        
        % V_gradient = x;
		%{
		x_1 = x( 1, 1 );
		x_2 = x( 1, 2 );
		% x_3 = x( 1, 3 );
		r = sqrt( x_1^2 + x_2^2 );
		%}
		x_1 = x( :, 1 );
		x_2 = x( :, 2 );
		r = sqrt( x_1.^2 + x_2.^2 );
		
		dV_dx1 = x_1 + gamma * ( x_2 .* cos( x_1 .* x_2 ) .* chi_r_func( r, R_minus_2, alpha_chi ) + sin( x_1 .* x_2 ) .* dchi_dr_func( r, R_minus_2, alpha_chi ) .* x_1 ./ r );
		dV_dx2 = alpha * x_2 + gamma * ( x_1 .* cos( x_1 .* x_2 ) .* chi_r_func( r, R_minus_2, alpha_chi ) + sin( x_1 .* x_2 ) .* dchi_dr_func( r, R_minus_2, alpha_chi ) .* x_2 ./ r );
		% dV_dx1 = x_1 + gamma * x_2 * cos( x_1 * x_2 );
		% dV_dx2 = alpha * x_2 + gamma * x_1 * cos( x_1 * x_2 );
		% dV_dx3 = x_3;
		
		% V_gradient = [ dV_dx1, dV_dx2, dV_dx3 ];
		V_gradient = [ dV_dx1, dV_dx2 ];
        
    end
    
    
    
    
    
    
    
    
    
    
    
    