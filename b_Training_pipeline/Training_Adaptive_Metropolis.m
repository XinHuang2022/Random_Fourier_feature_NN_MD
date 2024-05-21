 		
	function Training_status = Training_Adaptive_Metropolis( dim, N, lambda_1, lambda_2, lambda_3, C_const_bound, K_values, Num_parallel_worker, Num_replica, Training_data_store_path )	
		
		N_test = N;   					% N_test is the testing data set size

		M = 50;      					% M is the number of Adaptive Metropolis and Gradient Descent iterations
		delta_step = ( 2.4^2 ) / dim;			
		delta = delta_step / 10;		% delta is the step length in Metropolis trial proposals
		gamma_expo = 3 * dim - 2; 		% gamma_expo is the exponent parameter in Metropolis test.
		delta_lr_omega = 0.02;			% delta_lr is the learning rate for SGD with training on the frequencies
		
		num_K_values = length( K_values );
		Q = Num_replica;       					% Q is the number of replica for evaluating the error bar of the loss function

		Errors_K_rec_1 = zeros( num_K_values, 2 );    % Errors_K_rec_1 stores the mean and error bars of loss in the adaptive Metropolis algorithm
		Errors_K_rec_2 = zeros( num_K_values, 2 );    % Errors_K_rec_2 stores the mean and error bars of loss in the fixed random feature algorithm
		
		Error_rep_q_1 = zeros( num_K_values, Q );
		Error_rep_q_2 = zeros( num_K_values, Q );
		
		data_subfolder_name = fullfile( Training_data_store_path, 'training_data_temp' );
		
		early_stopping_steps_num = 20;
		parpool( round( Num_parallel_worker / 2 ) );
		tic 

		parfor q = 1 : Q
			
			% Load the q:th part of the dataset
			loaded_training_data_part_q = load( fullfile( data_subfolder_name, sprintf( 'Training_data_part_%d.mat', q ) ) ); % Load the corresponding part
			training_data_part_q = loaded_training_data_part_q.training_data_part;

			data_indices_for_training_q = 1 : N;
			x_data = training_data_part_q( data_indices_for_training_q, 1 : 2 );
			v_data = training_data_part_q( data_indices_for_training_q, 3 );
			v_prime_1_data = training_data_part_q( data_indices_for_training_q, 4 );
			v_prime_2_data = training_data_part_q( data_indices_for_training_q, 5 );
			y_data = [ v_data; v_prime_1_data; v_prime_2_data ];
			
			data_indices_for_testing_q = N + 1 : 2 * N;
			x_data_test = training_data_part_q( data_indices_for_testing_q, 1 : 2 );
			v_data_test = training_data_part_q( data_indices_for_testing_q, 3 );
			v_prime_1_data_test = training_data_part_q( data_indices_for_testing_q, 4 );
			v_prime_2_data_test = training_data_part_q( data_indices_for_testing_q, 5 );
			y_data_test = [ v_data_test; v_prime_1_data_test; v_prime_2_data_test ];
			
			N_num_Newton_beta = 50;
			
			for j = 1 : 1 : num_K_values
				
				K = K_values( 1, j );
				% rng(126)
				omega_AM_0 = 3 * randn( K, dim );
				omega_GD_0 = 3 * randn( K, dim );
				
				omega = omega_AM_0;
				S_v_mat = exp( 1i * ( x_data * omega' ) );
				omega_dim_1 = 1i * omega( :, 1 );
				S_v_grad_1_mat = ( omega_dim_1.' ) .* S_v_mat;
				omega_dim_2 = 1i * omega( :, 2 );
				S_v_grad_2_mat = ( omega_dim_2.' ) .* S_v_mat;
				S_mat = [ S_v_mat; S_v_grad_1_mat; S_v_grad_2_mat ];
				beta_hat_init = ( S_mat' * S_mat + lambda_1 * N * diag( ones( K, 1 ) ) ) \ ( S_mat' * y_data );
				beta_hat = Newton_Raphson_beta( S_mat, beta_hat_init, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );

				S_v_mat = [];  S_v_grad_1_mat = [];  S_v_grad_2_mat = [];  S_mat = [];
				
				Loss_min_rec_AM = 100;
				Loss_min_index_AM = 0;
				Loss_min_beta_rec_AM = ones( K, 1 );
				Loss_min_omega_rec_AM = ones( K, dim );
				
				for m = 1 : 1 : M
					
					r_normal_step = randn( K, dim );
					omega_prime = omega + delta * r_normal_step;
					S_v_mat_prime = exp( 1i * ( x_data * omega_prime' ) );
					omega_prime_dim_1 = 1i * omega_prime( :, 1 );
					S_v_grad_1_mat_prime = ( omega_prime_dim_1.' ) .* S_v_mat_prime;
					omega_prime_dim_2 = 1i * omega_prime( :, 2 );
					S_v_grad_2_mat_prime = ( omega_prime_dim_2.' ) .* S_v_mat_prime;
					S_mat_prime = [ S_v_mat_prime; S_v_grad_1_mat_prime; S_v_grad_2_mat_prime ];
					beta_hat_prime = ( S_mat_prime' * S_mat_prime + lambda_1 * N * diag( ones( K, 1 ) ) ) \ ( S_mat_prime' * y_data );
					beta_hat_prime = Newton_Raphson_beta( S_mat_prime, beta_hat_prime, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
					beta_hat_prime = Prox_Grad_beta( S_mat_prime, beta_hat_prime, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
					
					S_v_mat_prime = [];  S_v_grad_1_mat_prime = [];  S_v_grad_2_mat_prime = [];  S_mat_prime = [];
	  
					r_uniform = rand( K, 1 );
					beta_k_increase = r_uniform < ( ( beta_hat_prime ./ beta_hat ).^gamma_expo );
					omega_new = beta_k_increase .* omega_prime + ( 1 - beta_k_increase ) .* omega;
					S_v_mat_new = exp( 1i * ( x_data * omega_new' ) );
					omega_new_dim_1 = 1i * omega_new( :, 1 );
					S_v_grad_1_mat_new = ( omega_new_dim_1.' ) .* S_v_mat_new;
					omega_new_dim_2 = 1i * omega_new( :, 2 );
					S_v_grad_2_mat_new = ( omega_new_dim_2.' ) .* S_v_mat_new;
					S_mat_new = [ S_v_mat_new; S_v_grad_1_mat_new; S_v_grad_2_mat_new ];
					beta_hat_new = ( S_mat_new' * S_mat_new + lambda_1 * N * diag( ones( K, 1 ) ) ) \ ( S_mat_new' * y_data );
					beta_hat_new = Newton_Raphson_beta( S_mat_new, beta_hat_new, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
					beta_hat_new = Prox_Grad_beta( S_mat_new, beta_hat_new, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
					
					S_v_mat_new = [];  S_v_grad_1_mat_new =[];  S_v_grad_2_mat_new = [];  S_mat_new = [];
					
					S_v_mat_test = exp( 1i * ( x_data_test * omega_new' ) );
					S_v_grad_1_mat_test = ( omega_new_dim_1.' ) .* S_v_mat_test;
					S_v_grad_2_mat_test = ( omega_new_dim_2.' ) .* S_v_mat_test;
					S_mat_test = [ S_v_mat_test; S_v_grad_1_mat_test; S_v_grad_2_mat_test ];
					Loss_1 = ( 1 / N_test ) * norm( S_mat_test * beta_hat_new - y_data_test )^2;

					General_loss_AM_m = Loss_1;
					
					S_v_mat_test = [];  S_v_grad_1_mat_test = [];  S_v_grad_2_mat_test = [];  S_mat_test = [];
					
					if( General_loss_AM_m < Loss_min_rec_AM )
						Loss_min_index_AM = m;
						Loss_min_beta_rec_AM = beta_hat_new;
						Loss_min_omega_rec_AM = omega_new;
						Loss_min_rec_AM = General_loss_AM_m;
					end
					
					omega = omega_new;
					beta_hat = beta_hat_new;

					if( m > Loss_min_index_AM + early_stopping_steps_num )		% If we do not see a decrease in the generalized loss value for 5 steps of updates in \omega
						omega = Loss_min_omega_rec_AM;
						beta_hat = Loss_min_beta_rec_AM;
						break                          							% Then we implement the early-stopping
					end
					
				end
				
				omega_final = omega;
				beta_hat_final = beta_hat;
				
				omega_final_dim_1 = 1i * omega_final( :, 1 );
				omega_final_dim_2 = 1i * omega_final( :, 2 );
				
				S_v_mat_test = exp( 1i * ( x_data_test * omega_final' ) );
				S_v_grad_1_mat_test = ( omega_final_dim_1.' ) .* S_v_mat_test;
				S_v_grad_2_mat_test = ( omega_final_dim_2.' ) .* S_v_mat_test;
				S_mat_test = [ S_v_mat_test; S_v_grad_1_mat_test; S_v_grad_2_mat_test ];
				
				loss_K_adaptive = sum( abs( S_mat_test * beta_hat_final - y_data_test ).^2, 1 );
				Error_rep_q_1( j, q ) = loss_K_adaptive / N_test;
				
				S_v_mat_test = [];  S_v_grad_1_mat_test = [];  S_v_grad_2_mat_test = [];  S_mat_test = [];
				
				Freq_store_path = fullfile( 'b_Training_pipeline', 'frequencies' );
				if ~exist( Freq_store_path, 'dir' )			% Create 'frequencies' folder if it doesn't exist
					mkdir( Freq_store_path );
				end
				Ampl_store_path = fullfile( 'b_Training_pipeline', 'amplitudes' );
				if ~exist( Ampl_store_path, 'dir' )			% Create 'amplitudes' folder if it doesn't exist
					mkdir( Ampl_store_path );
				end
				
				filename_omega_AM = sprintf( 'frequency_omega_AM_K=%d_q=%d.csv', K, q );
				filepath_omega_AM = fullfile( Freq_store_path, filename_omega_AM );  % Full path including 'frequencies' subfolder
				writematrix( omega', filepath_omega_AM );
				beta_AM_real = real( beta_hat_final );
				beta_AM_imag = imag( beta_hat_final );
				beta_AM_final = [ beta_AM_real, beta_AM_imag ];
				filename_eta_AM = sprintf( 'amplitude_eta_AM_K=%d_q=%d.csv', K, q );
				filepath_eta_AM = fullfile( Ampl_store_path, filename_eta_AM );  % Full path including 'amplitudes' subfolder
				writematrix( beta_AM_final', filepath_eta_AM );
				
				
				
				omega_gd = omega_GD_0;
				S_v_mat_gd = exp( 1i * ( x_data * omega_gd' ) );
				omega_dim_1_gd = 1i * omega_gd( :, 1 );
				S_v_grad_1_mat_gd = ( omega_dim_1_gd.' ) .* S_v_mat_gd;
				omega_dim_2_gd = 1i * omega_gd( :, 2 );
				S_v_grad_2_mat_gd = ( omega_dim_2_gd.' ) .* S_v_mat_gd;
				S_mat_gd = [ S_v_mat_gd; S_v_grad_1_mat_gd; S_v_grad_2_mat_gd ];
				beta_hat_gd_init = ( S_mat_gd' * S_mat_gd + lambda_1 * N * diag( ones( K, 1 ) ) ) \ ( S_mat_gd' * y_data );
				beta_hat_gd_0 = Newton_Raphson_beta( S_mat_gd, beta_hat_gd_init, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
				beta_hat_gd = beta_hat_gd_0;

				Loss_1_init_gd = ( 1 / N ) * norm( S_mat_gd * beta_hat_gd - y_data )^2;
				Penalty_1_init_gd = lambda_1 * ( beta_hat_gd' * beta_hat_gd );
				Penalty_2_init_gd = lambda_2 * ( beta_hat_gd' * beta_hat_gd )^2;
				Penalty_3_init_gd = lambda_3 * max( sum( abs( real( beta_hat_gd ) ) + abs( imag( beta_hat_gd ) ) ) - C_const_bound, 0 );
				
				Loss_min_rec_GD = 100;
				Loss_min_index_GD = 0;

				for m_GD = 1 : 1 : M
					
					dL1_domega_1 = ( 2 / N ) * real( ( ( S_v_mat_gd * beta_hat_gd - v_data ).' * ( conj( S_v_mat_gd ) .* ( -1i * x_data( :, 1 ) ) ) ) .* ( beta_hat_gd.' ) );
					dL1_domega_2 = ( 2 / N ) * real( ( ( S_v_mat_gd * beta_hat_gd - v_data ).' * ( conj( S_v_mat_gd ) .* ( -1i * x_data( :, 2 ) ) ) ) .* ( beta_hat_gd.' ) );
					dL21_domega_1 = ( 2 / N ) * real( ( ( -1i ) * ( conj( beta_hat_gd ) .* S_v_mat_gd' ) + ( omega_gd( :, 1 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 1 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_1_mat_gd * beta_hat_gd - v_prime_1_data ) );
					dL21_domega_2 = ( 2 / N ) * real( ( ( omega_gd( :, 1 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 2 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_1_mat_gd * beta_hat_gd - v_prime_1_data ) );
					dL22_domega_1 = ( 2 / N ) * real( ( ( omega_gd( :, 2 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 1 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_2_mat_gd * beta_hat_gd - v_prime_2_data ) );
					dL22_domega_2 = ( 2 / N ) * real( ( ( -1i ) * ( conj( beta_hat_gd ) .* S_v_mat_gd' ) + ( omega_gd( :, 2 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 2 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_2_mat_gd * beta_hat_gd - v_prime_2_data ) );
					
					dL_domega_1 = dL1_domega_1.' + dL21_domega_1 + dL22_domega_1;
					dL_domega_2 = dL1_domega_2.' + dL21_domega_2 + dL22_domega_2;
					
					dL_domega_part_1 = [ dL_domega_1, dL_domega_2 ];
					dL_domega = dL_domega_part_1;
					
					omega_gd_new = omega_gd - delta_lr_omega * dL_domega;
					
					omega_gd = omega_gd_new;

					S_v_mat_gd = exp( 1i * ( x_data * omega_gd' ) );
					omega_dim_1_gd = 1i * omega_gd( :, 1 );
					S_v_grad_1_mat_gd = ( omega_dim_1_gd.' ) .* S_v_mat_gd;
					omega_dim_2_gd = 1i * omega_gd( :, 2 );
					S_v_grad_2_mat_gd = ( omega_dim_2_gd.' ) .* S_v_mat_gd;
					S_mat_gd = [ S_v_mat_gd; S_v_grad_1_mat_gd; S_v_grad_2_mat_gd ];
					beta_hat_gd = ( S_mat_gd' * S_mat_gd + lambda_1 * N * diag( ones( K, 1 ) ) ) \ ( S_mat_gd' * y_data );
					beta_hat_gd = Newton_Raphson_beta( S_mat_gd, beta_hat_gd, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
					beta_hat_gd = Prox_Grad_beta( S_mat_gd, beta_hat_gd, y_data, N, K, lambda_1, lambda_2, lambda_3, C_const_bound, N_num_Newton_beta );
					beta_hat_gd_0 = beta_hat_gd;
					
					S_v_mat_gd_test = exp( 1i * ( x_data_test * omega_gd' ) );
					S_v_grad_1_mat_gd_test = ( omega_dim_1_gd.' ) .* S_v_mat_gd_test;
					S_v_grad_2_mat_gd_test = ( omega_dim_2_gd.' ) .* S_v_mat_gd_test;
					S_mat_gd_test = [ S_v_mat_gd_test; S_v_grad_1_mat_gd_test; S_v_grad_2_mat_gd_test ];
					Loss_1_GD = ( 1 / N_test ) * norm( S_mat_gd_test * beta_hat_gd - y_data_test )^2;
					Penalty_1_GD = lambda_1 * ( beta_hat_gd' * beta_hat_gd );
					Penalty_2_GD = lambda_2 * ( beta_hat_gd' * beta_hat_gd )^2;
					Penalty_3_GD = lambda_3 * max( sum( abs( real( beta_hat_gd ) ) + abs( imag( beta_hat_gd ) ) ) - C_const_bound, 0 );
					General_loss_GD_m = Loss_1_GD + Penalty_1_GD + Penalty_2_GD + Penalty_3_GD;
					
					% clear S_v_mat_gd_test  S_v_grad_1_mat_gd_test  S_v_grad_2_mat_gd_test  S_mat_gd_test
					S_v_mat_gd_test = [];  S_v_grad_1_mat_gd_test = [];  S_v_grad_2_mat_gd_test = [];  S_mat_gd_test = [];
				   
					if( General_loss_GD_m < Loss_min_rec_GD )
						Loss_min_index_GD = m_GD;
						Loss_min_rec_GD = General_loss_GD_m;
					end

					if( m_GD > Loss_min_index_GD + early_stopping_steps_num )    	% If we do not see a decrease in the generalized loss value for 5 steps of updates in \omega
						break                          								% Then we implement the early-stopping
					end
					
				end

				S_v_mat_gd_test = exp( 1i * ( x_data_test * omega_gd' ) );
				S_v_grad_1_mat_gd_test = ( omega_dim_1_gd.' ) .* S_v_mat_gd_test;
				S_v_grad_2_mat_gd_test = ( omega_dim_2_gd.' ) .* S_v_mat_gd_test;
				S_mat_gd_test = [ S_v_mat_gd_test; S_v_grad_1_mat_gd_test; S_v_grad_2_mat_gd_test ];

				loss_K_gd = sum( abs( S_mat_gd_test * beta_hat_gd - y_data_test ).^2, 1 );
				Error_rep_q_2( j, q ) = loss_K_gd / N_test;
				
				% clear S_v_mat_gd_test  S_v_grad_1_mat_gd_test  S_v_grad_2_mat_gd_test  S_mat_gd_test
				S_v_mat_gd_test = [];  S_v_grad_1_mat_gd_test = [];  S_v_grad_2_mat_gd_test = [];  S_mat_gd_test = [];
				
				filename_omega_GD = sprintf( 'frequency_omega_GD_K=%d_q=%d.csv', K, q );
				filepath_omega_GD = fullfile( Freq_store_path, filename_omega_GD );  % Full path including 'frequencies' subfolder
				writematrix( omega_gd', filepath_omega_GD );
				beta_GD_real = real( beta_hat_gd );
				beta_GD_imag = imag( beta_hat_gd );
				beta_GD_final = [ beta_GD_real, beta_GD_imag ];
				filename_eta_GD = sprintf( 'amplitude_eta_GD_K=%d_q=%d.csv', K, q );
				filepath_eta_GD = fullfile( Ampl_store_path, filename_eta_GD );  % Full path including 'amplitudes' subfolder
				writematrix( beta_GD_final', filepath_eta_GD );

			end

			fprintf( 'Independent replica Number %d for training completed.\n', q );
			
		end
		poolobj = gcp( 'nocreate' );
		delete( poolobj );

		toc
		
		Errors_K_rec_1( :, 1 ) = mean( Error_rep_q_1, 2 );
		Errors_K_rec_1( :, 2 ) = std( Error_rep_q_1, 0, 2 );
		
		Errors_K_rec_2( :, 1 ) = mean( Error_rep_q_2, 2 );
		Errors_K_rec_2( :, 2 ) = std( Error_rep_q_2, 0, 2 );

		Loss_training_subfolder_name = fullfile( 'b_Training_pipeline', 'training_loss_record_temp' );
		% Create the subfolder if it doesn't exist
		if ~exist( Loss_training_subfolder_name, 'dir' )
			mkdir( Loss_training_subfolder_name );
		end
		% Save the current part to the subfolder
		save( fullfile( Loss_training_subfolder_name, sprintf( 'Training_loss_J=%d_Kmax=%d_R=%d_AM.mat', N, max( K_values ), Q ) ), 'Error_rep_q_1' );
		save( fullfile( Loss_training_subfolder_name, sprintf( 'Training_loss_J=%d_Kmax=%d_R=%d_GD.mat', N, max( K_values ), Q ) ), 'Error_rep_q_2' );
		save( fullfile( Loss_training_subfolder_name, sprintf( 'Training_loss_stat_J=%d_Kmax=%d_R=%d_AM.mat', N, max( K_values ), Q ) ), 'Errors_K_rec_1' );
		save( fullfile( Loss_training_subfolder_name, sprintf( 'Training_loss_stat_J=%d_Kmax=%d_R=%d_GD.mat', N, max( K_values ), Q ) ), 'Errors_K_rec_2' );
		
		% Previous data recorded with TensorFlow training for comparison
		% K_values_sgd = [ 64, 128, 256, 512, 1024 ];
		% Errors_K_rec_3 = [ 0.4995, 0.1517, 0.0910, 0.0556, 0.0309 ];
		
		fig3 = figure( 3 );
		errorbar( K_values, Errors_K_rec_1( :, 1 ), Errors_K_rec_1( :, 2 ) * 1.96 / sqrt( Q ), '-*' );
		hold on
		errorbar( K_values, Errors_K_rec_2( :, 1 ), Errors_K_rec_2( :, 2 ) * 1.96 / sqrt( Q ), '-o' );
		loglog( K_values, 25 * K_values.^(-1) );
		hold off
		legend( 'loss for adaptive Metropolis algorithm', 'loss for Gradient Descent on frequency', 'reference $\mathcal{O}(K^{-1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
		xlabel( '$K$', 'fontsize', 20, 'interpreter', 'latex' );
		ylabel( 'Loss on the test set', 'fontsize', 20, 'interpreter', 'latex' );
		title_string = sprintf( 'Log-log plot of testing loss, J=%d, d=%d, Q=%d', N, dim, Q );
		title( title_string, 'fontsize', 20 );
		% title( title_string, 'fontsize', 20, 'interpreter', 'latex' )
		set( gca, 'YScale', 'log');
		set( gca, 'XScale', 'log');
		
		Fig_FolderName = fullfile( 'b_Training_pipeline', 'Figures_save_training_loss' );
		if ~exist( Fig_FolderName, 'dir' )			% Create 'frequencies' folder if it doesn't exist
			mkdir( Fig_FolderName );
		end
		training_loss_figure_name = sprintf( 'Training_loss_plot_K=%d_J=%d_N_rep=%d.fig', max( K_values ), N, Q );
		saveas( fig3, fullfile( Fig_FolderName, training_loss_figure_name ) ); 
		
		%{
		figure( 2 )
		alpha = sqrt( 2 );
		gamma = 2;
		x1_plot_test = ( -5 : 0.01 : 5 )';
		x2_plot_test = x1_plot_test;
		x_plot_test = [ x1_plot_test, x2_plot_test ];
		r_plot_points = sqrt( sum( x_plot_test.^2, 2 ) );
		num_x_plot_point = length( x1_plot_test );
		v1_x_plot = zeros( num_x_plot_point, 1 );
		for i = 1 : 1 : num_x_plot_point
			r_point_i = r_plot_points( i, 1 );
			x1_i = x1_plot_test( i, 1 );
			x2_i = x2_plot_test( i, 1 );
			% v1_x_plot( i, 1 ) = ( x1_i^2 / 2 + alpha * x2_i^2 / 2 + gamma * sin( x1_i * x2_i ) * chi_r_func( r_point_i, R_minus_2, R_minus_1 ) ) * chi_r_func( r_point_i, R_0, R_1 );
			v1_x_plot( i, 1 ) = ( x1_i^2 / 2 + alpha * x2_i^2 / 2 + gamma * sin( x1_i * x2_i ) * chi_r_func( r_point_i, R_minus_2, alpha_chi ) ) * chi_r_func( r_point_i, R_c, alpha_chi );
		end
		
		S_mat_1_plot = exp( 1i * ( x_plot_test * omega' ) );
		alpha_1_output = S_mat_1_plot * beta_hat_final;
		
		S_mat_2_plot = exp( 1i * ( x_plot_test * omega_gd' ) );
		alpha_2_output = S_mat_2_plot * beta_hat_gd;
		
		plot( x1_plot_test, alpha_1_output );
		hold on
		plot( x1_plot_test, alpha_2_output );
		plot( x1_plot_test, v1_x_plot );
		%}
		%{
		A1_plot = 0 : 0.1 : 4 * C_const_bound;
		chi_plot = zeros( length( A1_plot ), 1 );
		dchi_dA_plot = zeros( length( A1_plot ), 1 );
		d2chi_dA2_plot = zeros( length( A1_plot ), 1 );
		for i_c = 1 : 1 : length( A1_plot )
			A1_ic = A1_plot( 1, i_c );
			chi_plot( i_c, 1 ) = chi_func( A1_ic, C_const_bound, epsilon_chi );
			dchi_dA_plot( i_c, 1 ) = dchi_dA_func( A1_ic, C_const_bound, epsilon_chi );
			d2chi_dA2_plot( i_c, 1 ) = d2chi_dA2_func( A1_ic, C_const_bound, epsilon_chi );
		end
		hold on
		plot( A1_plot, chi_plot );
		plot( A1_plot, dchi_dA_plot );
		plot( A1_plot, d2chi_dA2_plot );
		legend( '\chi(A)', 'd\chi/dA', 'd^2\chi/dA^2' )
		%}
		
		%{
		figure( 2 )
		plot( 1 : m, General_loss_rec_update_omega( 1 : m, 1 ), 'o-' );
		xlabel( 'Number of updates in \omega' );
		ylabel( 'Generalized loss value' );
		title( 'Adaptive Metropolis method for updates in \omega' )
		
		
		figure( 3 )
		plot( 1 : m_GD, General_loss_rec_update_omega( 1 : m_GD, 2 ), '*-' );
		xlabel( 'Number of updates in \omega' );
		ylabel( 'Generalized loss value' );
		title( 'Gradient Descent method for updates in \omega' )
		%}
		
		Training_status = 1;
		
	end	
