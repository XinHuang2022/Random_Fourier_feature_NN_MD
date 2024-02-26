	
	function beta_hat_new = Newton_Raphson_beta( S_mat, beta_hat_0, y_data, N, K, lambda_1, lambda_2, N_num_Newton_beta )
	
		beta_vec = beta_hat_0;
		error_rec = zeros( N_num_Newton_beta, 1 );
		delta_beta = 1;
		Tol = 10^(-8);
        for m = 1 : 1 : N_num_Newton_beta
			
			% J_beta = ( S_mat' * S_mat ) + 1 * lambda_1 * N * diag( ones( K, 1 ) );
			% F_beta = ( ( S_mat' * S_mat ) + 1 * lambda_1 * N * diag( ones( K, 1 ) ) ) * beta_vec - S_mat' * y_data;
            F_beta = ( S_mat' * S_mat + ( 2 * lambda_1 * N ) * diag( ones( K, 1 ) ) ) * beta_vec + 4 * lambda_2 * N * norm( beta_vec )^2 * beta_vec - S_mat' * y_data;
			J_beta = S_mat' * S_mat + ( 2 * lambda_1 * N + 4 * lambda_2 * N * ( beta_vec' * beta_vec ) ) * diag( ones( K, 1 ) ) + 8 * lambda_2 * N * ( beta_vec * beta_vec.' );
			delta_beta = -J_beta \ F_beta;
            error_rec( m, 1 ) = norm( delta_beta );
			beta_vec = beta_vec + delta_beta;
			
			if( norm( delta_beta ) < Tol )
				break;
			end
			
		end
		beta_hat_new = beta_vec;
		
		%{
		figure( 2 )
		% plot( 1 : N_num_Newton_beta, error_rec, '.-' );
		% semilogy( 1 : N_num_Newton_beta, error_rec, '.-' );
		semilogy( 1 : m, error_rec( 1 : m, 1 ), '.-' );
		%}
		
	end
	
	
	
	
	