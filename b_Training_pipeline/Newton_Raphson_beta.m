	
	function beta_hat_new = Newton_Raphson_beta( S_mat, beta_hat_0, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const, N_num_Newton_beta )
	
		beta_vec = beta_hat_0;
		% error_rec = zeros( N_num_Newton_beta, 1 );
		Abs_Tol = 10^(-6);
        Rel_Tol = 1 * 10^(-3);

		Loss_1_init = ( 1 / N ) * ( S_mat * beta_hat_0 - y_data )' * ( S_mat * beta_hat_0 - y_data );
        Penalty_1_init = lambda_1 * ( beta_hat_0' * beta_hat_0 );
        Penalty_2_init = lambda_2 * ( beta_hat_0' * beta_hat_0 )^2;
		Penalty_3_init = lambda_3 * max( sum( abs( real( beta_hat_0 ) ) + abs( real( beta_hat_0 ) ) ) - C_const, 0 );
        generalize_loss_init = Loss_1_init + Penalty_1_init + Penalty_2_init + Penalty_3_init;
		
		S_design_mat = S_mat' * S_mat;
        generalize_loss_prev = generalize_loss_init;

        for m = 1 : 1 : N_num_Newton_beta
			
			% F_beta is the gradient of the loss function w.r.t. the amplitude parameters
			% J_beta is the Jacobian of F_beta, i.e., it is the Hessian of the loss function w.r.t. the amplitude parameters
			
			a_vec = real( beta_vec );
			b_vec = imag( beta_vec );
			
			beta_vec_L2_square = beta_vec' * beta_vec;
			F_a_vec_1 = ( 2 / N ) * real( S_design_mat * beta_vec ) - ( 2 / N ) * real( S_mat' * y_data ) + 2 * lambda_1 * a_vec + 4 * lambda_2 * beta_vec_L2_square * a_vec;
			F_b_vec_1 = ( 2 / N ) * imag( S_design_mat * beta_vec ) - ( 2 / N ) * imag( S_mat' * y_data ) + 2 * lambda_1 * b_vec + 4 * lambda_2 * beta_vec_L2_square * b_vec;
			
			J_a_mat_1 = ( 2 / N ) * real( S_design_mat ) + ( 2 * lambda_1 + 4 * lambda_2 * beta_vec_L2_square ) * eye( K ) + 8 * lambda_2 * ( a_vec * a_vec' );
			J_b_mat_1 = ( 2 / N ) * real( S_design_mat ) + ( 2 * lambda_1 + 4 * lambda_2 * beta_vec_L2_square ) * eye( K ) + 8 * lambda_2 * ( b_vec * b_vec' );
			
            % F_a_vec_1( 1, 1 )
            % J_b_mat_1( 2, 1 )
            %{
			F_a_vec_2 = alpha_3 * ( a_vec ./ abs( beta_vec ) );
			F_b_vec_2 = alpha_3 * ( b_vec ./ abs( beta_vec ) );
			
			J_a_mat_2 = alpha_3 * diag( ( b_vec.^2 ) ./ ( abs( beta_vec ).^3 ) );
			J_b_mat_2 = alpha_3 * diag( ( a_vec.^2 ) ./ ( abs( beta_vec ).^3 ) );
            

			F_a_vec = F_a_vec_1 + F_a_vec_2;
			J_a_mat = J_a_mat_1 + J_a_mat_2;
			F_b_vec = F_b_vec_1 + F_b_vec_2;
			J_b_mat = J_b_mat_1 + J_b_mat_2;
            %}
			delta_a = -J_a_mat_1 \ F_a_vec_1;
			delta_b = -J_b_mat_1 \ F_b_vec_1;
            % delta_a = -J_a_mat_2 \ F_a_vec_2;
			% delta_b = -J_b_mat_2 \ F_b_vec_2;
			% delta_a = -J_a_mat \ F_a_vec;
			% delta_b = -J_b_mat \ F_b_vec;
			delta_beta = delta_a + 1i * delta_b;
		
            beta_vec_new = beta_vec + delta_beta;
            % error_rec( m, 1 ) = sqrt( norm( delta_a )^2 + norm( delta_b )^2 );
			% error_rec( m, 1 ) = norm( delta_beta );
            Loss_1 = ( 1 / N ) * norm( S_mat * beta_vec_new - y_data )^2;
            Penalty_1 = lambda_1 * ( beta_vec_new' * beta_vec_new );
            Penalty_2 = lambda_2 * ( beta_vec_new' * beta_vec_new )^2;
			Penalty_3 = lambda_3 * max( sum( abs( real( beta_vec_new ) ) + abs( real( beta_vec_new ) ) ) - C_const, 0 );
            stop_check =  abs( Loss_1 + Penalty_1 + Penalty_2 + Penalty_3 - generalize_loss_prev );
			
            beta_vec = beta_vec_new;
            
			%  if( norm( delta_beta ) / sqrt( K ) < Tol )
            
            if( ( stop_check < Abs_Tol ) || ( stop_check / generalize_loss_prev < Rel_Tol ) )
                % fprintf( 'Newton-Raphson method converged with %d iterations\n', m );
				break;
            end

            if( m == N_num_Newton_beta )
                fprintf( 'Newton-Raphson method uses %d iterations\n', m );
            end
            
            generalize_loss_prev = Loss_1 + Penalty_1 + Penalty_2 + Penalty_3;

		end
		beta_hat_new = beta_vec;
		% error_rec
		
		%{
		figure( 2 )
		% plot( 1 : N_num_Newton_beta, error_rec, '.-' );
		% semilogy( 1 : N_num_Newton_beta, error_rec, '.-' );
		semilogy( 1 : m, error_rec( 1 : m, 1 ), '.-' );
		%}
		
	end
	
	
	
	
	