	
	function beta_hat_new = Prox_Grad_beta( S_mat, beta_hat_0, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const, N_num_Newton_beta )
	
		beta_vec = beta_hat_0;
		Rel_Tol = 1 * 10^(-3);

        Loss_1_init = ( 1 / N ) * ( S_mat * beta_hat_0 - y_data )' * ( S_mat * beta_hat_0 - y_data );
        Penalty_1_init = lambda_1 * ( beta_hat_0' * beta_hat_0 );
        Penalty_2_init = lambda_2 * ( beta_hat_0' * beta_hat_0 )^2;
		Penalty_3_init = lambda_3 * max( sum( abs( real( beta_hat_0 ) ) + abs( real( beta_hat_0 ) ) ) - C_const, 0 );
        generalize_loss_init = Loss_1_init + Penalty_1_init + Penalty_2_init + Penalty_3_init;
		% Loss_rec = zeros( N_num_Newton_beta, 1 );
		% Loss_rec( 1, 1 ) = generalize_loss_init;
		
		S_design_mat = S_mat' * S_mat;
		gamma = 0.5 / K;
		beta_decrease_factor = 0.95;

        a_vec = real( beta_vec );
		b_vec = imag( beta_vec );
        generalize_loss_prev = generalize_loss_init;
		beta_vec_old = beta_vec;
		
        for m = 1 : 1 : N_num_Newton_beta
			
			% F_beta is the gradient of the loss function w.r.t. the amplitude parameters
			% J_beta is the Jacobian of F_beta, i.e., it is the Hessian of the loss function w.r.t. the amplitude parameters
			
			sum_ampli_L1 = sum( abs( a_vec ) + abs( b_vec ) );
			if( sum_ampli_L1 > C_const )
			    
				%{
                if( mod( m, 10 ) == 0 )
                    fprintf( 'K = %d, iteration m = %d, p3-term > 0, using subgradient method\n', K, m );
                end
				%}
				
				extrapol_para_m = m / ( m + 3 );
				eta_vec = beta_vec + extrapol_para_m * ( beta_vec - beta_vec_old );
				real_eta_vec = real( eta_vec );
				imag_eta_vec = imag( eta_vec );
				
				eta_vec_L2_square = eta_vec' * eta_vec;
				gradient_Loss_1_eta = ( 2 / N ) * ( S_design_mat * eta_vec - S_mat' * y_data );
				
				F_grad_eta_real = real( gradient_Loss_1_eta ) + 2 * lambda_1 * real_eta_vec + 4 * lambda_2 * eta_vec_L2_square * real_eta_vec;
				w_input_vec_eta_real = real_eta_vec - gamma * F_grad_eta_real;
				a_vec_new = Soft_threshold_func( w_input_vec_eta_real, lambda_3 * gamma, K );
				
				F_grad_eta_imag = imag( gradient_Loss_1_eta ) + 2 * lambda_1 * imag_eta_vec + 4 * lambda_2 * eta_vec_L2_square * imag_eta_vec;
				w_input_vec_eta_imag = imag_eta_vec - gamma * F_grad_eta_imag;
				b_vec_new = Soft_threshold_func( w_input_vec_eta_imag, lambda_3 * gamma, K );
				
				a_vec = a_vec_new;
				b_vec = b_vec_new;
				beta_vec_new = a_vec + 1i * b_vec;
			
			else
			    
			    beta_vec_L2_square = beta_vec' * beta_vec;
				F_a_vec_1 = ( 2 / N ) * real( S_design_mat * beta_vec ) - ( 2 / N ) * real( S_mat' * y_data ) + 2 * lambda_1 * a_vec + 4 * lambda_2 * beta_vec_L2_square * a_vec;
				F_b_vec_1 = ( 2 / N ) * imag( S_design_mat * beta_vec ) - ( 2 / N ) * imag( S_mat' * y_data ) + 2 * lambda_1 * b_vec + 4 * lambda_2 * beta_vec_L2_square * b_vec;
				
				J_a_mat_1 = ( 2 / N ) * real( S_design_mat ) + ( 2 * lambda_1 + 4 * lambda_2 * beta_vec_L2_square ) * eye( K ) + 8 * lambda_2 * ( a_vec * a_vec' );
				J_b_mat_1 = ( 2 / N ) * real( S_design_mat ) + ( 2 * lambda_1 + 4 * lambda_2 * beta_vec_L2_square ) * eye( K ) + 8 * lambda_2 * ( b_vec * b_vec' );

				delta_a = -J_a_mat_1 \ F_a_vec_1;
			    delta_b = -J_b_mat_1 \ F_b_vec_1;
				a_vec = a_vec + delta_a;
				b_vec = b_vec + delta_b;
				beta_vec_new = a_vec + 1i * b_vec;
			
			end
  
            Loss_1 = ( 1 / N ) * norm( S_mat * beta_vec_new - y_data )^2;
            Penalty_1 = lambda_1 * ( beta_vec_new' * beta_vec_new );
            Penalty_2 = lambda_2 * ( beta_vec_new' * beta_vec_new )^2;
			Penalty_3 = lambda_3 * max( sum( abs( real( beta_vec_new ) ) + abs( imag( beta_vec_new ) ) ) - C_const, 0 );
            stop_check =  abs( Loss_1 + Penalty_1 + Penalty_2 + Penalty_3 - generalize_loss_prev ) / generalize_loss_prev;
            % if( mod( m, 10 ) == 0 )
                % fprintf( 'K = %d, Iteration %d, Loss_1 = %4.3e, Penalty_1 = %4.3e, Penalty_2 = %4.3e, Penalty_3 = %4.3e, stop_check = %4.3e \n', K, m, Loss_1, Penalty_1, Penalty_2, Penalty_3, stop_check );
            % end

			beta_vec_old = beta_vec;
			beta_vec = beta_vec_new;
            
			% if( sqrt( sum( delta_a.^2 + delta_b.^2 ) )  < Tol )
            if( stop_check < Rel_Tol )
                % fprintf( 'Fixed-point method converged with %d iterations\n', m );
				break;
            end
            %{
            if( m == N_num_Newton_beta )
                fprintf( 'Fixed-point method uses %d iterations\n', m );
            end
            %}
            generalize_loss_prev = Loss_1 + Penalty_1 + Penalty_2 + Penalty_3;
			% Loss_rec( m, 1 ) = generalize_loss_prev;
			gamma = gamma * beta_decrease_factor;

		end

		beta_hat_new = beta_vec;

		
		%{
		figure( 2 )
		% plot( 1 : N_num_Newton_beta, error_rec, '.-' );
		% semilogy( 1 : N_num_Newton_beta, error_rec, '.-' );
		semilogy( 1 : m, error_rec( 1 : m, 1 ), '.-' );
		% plot( 1 : m, error_rec( 1 : m, 1 ), '.-' );
        %}
		%{
        if( K == 2048 )
		    figure( 3 )
		    % plot( 1 : N_num_Newton_beta, error_rec, '.-' );
		    % semilogy( 1 : N_num_Newton_beta, error_rec, '.-' );
		    semilogy( 1 : m, Loss_rec( 1 : m, 1 ), '.-' );
        end
		%}
		
	end
	
	
	
	
	