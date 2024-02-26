    
	function V_approx_2d_input = V_approx_2d_plot( x1, x2, R_a, R_b, sigma_a, sigma_b )
		
		omega = readmatrix( './frequencies/frequency_omega_AM_K=512_q=32.csv' );
		% eta = readmatrix( './parameter_set/K=256_paraset/weight_parameters_K=256.txt' );
		% eta = readmatrix( './parameter_set/K=1024_paraset/weight_parameters_K=1024_dt_2e-4.txt' );
		eta = readmatrix( './amplitudes/amplitude_eta_AM_K=512_q=32.csv' );
		eta_re = eta( 1, : );
		eta_im = eta( 2, : );
		
		V_int_x_i = zeros( 'like', x1 );
        for k = 1 : 1 : length( eta_re )
            omega_k_1 = omega( 1, k );
			omega_k_2 = omega( 2, k );
			inner_prod = x1 * omega_k_1 + x2 * omega_k_2;
            V_int_summand = eta_re( 1, k ) * cos( inner_prod ) - eta_im( 1, k ) * sin( inner_prod );
            V_int_x_i = V_int_x_i + V_int_summand;
		end
		
		alpha = sqrt( 2 );
		gamma = 2;
		
		%{
		A_mat = [ 3 * R_b^2, 2 * R_b, 1, 0; 3 * R_a^2, 2 * R_a, 1, 0; R_b^3, R_b^2, R_b, 1; R_a^3, R_a^2, R_a, 1 ];
		b_vec = [ 0; 0; 0; 1 ];
		c_para = A_mat \ b_vec;
		a = c_para( 1, 1 );  b = c_para( 2, 1 );  c = c_para( 3, 1 );  d = c_para( 4, 1 );
		%}
		
		R = sqrt( x1.^2 + x2.^2 );
		
        % R_minus_1 = 2;
		% V_true_outside_R0 = x1.^2 / 2 + alpha * x2.^2 / 2;
		V_true_outside_R0 = x1.^2 / 2 + alpha * x2.^2 / 2 + gamma * sin( x1 .* x2 ) .* exp( -( R - R_a ).^2 / ( 2 * sigma_a ) );
		
		% alpha_chi_2 = 0.1 * alpha_chi;
		% V_approx_2d_input = ( R < R_a ) .* V_int_x_i ...
							% + ( ( R >= R_a ) .* ( R <= R_b ) ) .* ( V_int_x_i .* exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) + V_true_outside_R0 .* ( 1 - exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) ) )...
							% + ( R > R_b ) .* V_true_outside_R0;
         % + ( ( R >= R_a ) ) .* ( V_int_x_i .* exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) + V_true_outside_R0 .* ( 1 - exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) ) );
                            % + ( ( R >= R_a ) .* ( R <= R_b ) ) .* ( V_int_x_i .* exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) + V_true_outside_R0 .* ( 1 - exp( -( R - R_a ).^2 / ( 2 * alpha_chi ) ) ) )...
		V_approx_2d_input = ( R < R_b ) .* V_int_x_i ...
							+ ( R >= R_b ) .* ( V_int_x_i .* exp( -( R - R_b ).^2 / ( 2 * sigma_b ) ) + V_true_outside_R0 .* ( 1 - exp( -( R - R_b ).^2 / ( 2 * sigma_b ) ) ) );
							% + ( R > R_b ) .* V_true_outside_R0;					
							
		% V_approx_2d_input = V_int_x_i;
    end
		
		
		
	
	
	
	
	
	