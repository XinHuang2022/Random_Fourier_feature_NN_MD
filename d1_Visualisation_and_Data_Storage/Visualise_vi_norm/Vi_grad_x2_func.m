    function dVi_dx2 = Vi_grad_x2_func( x1, x2, alpha, gamma, R_a, sigma_a, R_b, sigma_b )
		
		R = sqrt( x1.^2 + x2.^2 );
		dVi_dx2 = chi_r_func( R, R_b, sigma_b ) .* dv_dx2_func( x1, x2, alpha, gamma, R_a, sigma_a )...
				  + V_true_2d_func( x1, x2, alpha, gamma, R_a, sigma_a ) .* dchi_dr_func( R, R_b, sigma_b ) .* x2 ./ R;
	
	end