    
	function dVi_dx1 = Vi_grad_x1_func( x1, x2, alpha, gamma, R_a, sigma_a, R_b, sigma_b )
		
		R = sqrt( x1.^2 + x2.^2 );
		dVi_dx1 = chi_r_func( R, R_b, sigma_b ) .* dv_dx1_func( x1, x2, alpha, gamma, R_a, sigma_a )...
				  + V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a ) .* dchi_dr_func( R, R_b, sigma_b ) .* x1 ./ R;
	
	end
		
		
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	