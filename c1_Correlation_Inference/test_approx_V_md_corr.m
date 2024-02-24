    
	function corr_curve_info = test_approx_V_md_corr( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, Num_parallel_worker )
	
		%{
		% First load the model trained with Tensorflow and Keras
		modelFolder_2 = './Model_save_HO_shallow_1';
		net = importTensorFlowNetwork( modelFolder_2, 'TargetNetwork', 'dlnetwork' );
		analyzeNetwork(net)
		%}
		
		% Load the reference correlation curve obtained from the true potential
		% load( 'correlation_md_curve_h_0.0005_Msp=2_18.mat' );
		load( 'correlation_md_curve_h_0.00500_beta=1_Msp=2_28.mat' );
		Num_tau_points = round( tau / Delta_tau_plot ) + 1;
		correlation_curve_ref_x = corr_curve_info_1( 1, 1 : Num_tau_points );
		correlation_curve_ref_p = corr_curve_info_1( 4, 1 : Num_tau_points );
		
		beta = 1;
		dim = 2;
		% x_corr = 1;
		alpha = sqrt( 2 );
		gamma = 2;
		
		R_a = 2;    
		R_b = 3;    % R_0 = 4;
		sigma_a = 1;
		sigma_b = 0.2;
		
		% K_values = [ 8, 16, 32, 64, 128, 256, 512, 1024 ];
		K_values = [ 8, 16, 32, 64, 128, 256, 512 ];
		% K_values = [ 8, 16, 32, 64, 128, 256 ];
		% K_values = 64;
		num_K_value = length( K_values );
		
		L1_diff_rec_AM_x = zeros( num_K_value, 2 );    % first row for storing the mean of the L1-diff with R replica, second row for the std of the L1-diff with R replica.
		L1_diff_rec_GD_x = zeros( num_K_value, 2 );
		L1_diff_rec_AM_p = zeros( num_K_value, 2 );
		L1_diff_rec_GD_p = zeros( num_K_value, 2 );

		% Num_parallel_worker = 128;
		% Num_parallel_worker = 8;
		parpool( Num_parallel_worker );
		L_x_vectorize = 1;
		
		
		for i_a = 1 : 1 : num_K_value
		
			K = K_values( 1, i_a );

			% tau = 2;
			% Delta_tau_plot = 0.1;
			
			% M_sp = 10^5;
			% M_sp = 2^18;
			% M_sp = 2^12;
			Correlation_samples_AM_x = zeros( M_sp, Num_tau_points );
			Correlation_samples_AM_p = zeros( M_sp, Num_tau_points );
			Correlation_samples_GD_x = zeros( M_sp, Num_tau_points );
			Correlation_samples_GD_p = zeros( M_sp, Num_tau_points );
			
			% x_0_sample_store_AM = zeros( M_sp, dim );
			% x_0_sample_store_GD = zeros( M_sp, dim );
			h_x_equi = 0.1;
			% h_x_sample = 0.0005;
			
			% profile on
			% Num_parallel_worker = 10;
			
			M_sp_on_worker = M_sp / Num_parallel_worker;
			
			correlation_md_record_AM_x = zeros( Num_replica, Num_tau_points );
			correlation_md_record_AM_p = zeros( Num_replica, Num_tau_points );
			correlation_md_record_GD_x = zeros( Num_replica, Num_tau_points );
			correlation_md_record_GD_p = zeros( Num_replica, Num_tau_points );
			
			tic
			for R = 1 : Num_replica
				% seed_a = R^2;
				% rng( seed_a );
				
				filename_omega_GD = sprintf( 'frequency_omega_GD_K=%d_q=%d.csv', K, R );
				filepath_omega_GD = fullfile( 'frequencies', filename_omega_GD );  % Full path including 'frequencies' subfolder
				filename_omega_AM = sprintf( 'frequency_omega_AM_K=%d_q=%d.csv', K, R );
				filepath_omega_AM = fullfile( 'frequencies', filename_omega_AM );  % Full path including 'frequencies' subfolder
				
				filename_eta_GD = sprintf( 'amplitude_eta_GD_K=%d_q=%d.csv', K, R );
				filepath_eta_GD = fullfile( 'amplitudes', filename_eta_GD );  % Full path including 'amplitudes' subfolder
				filename_eta_AM = sprintf( 'amplitude_eta_AM_K=%d_q=%d.csv', K, R );
				filepath_eta_AM = fullfile( 'amplitudes', filename_eta_AM );  % Full path including 'amplitudes' subfolder
				
				% First load the model parameters
				% omega = readmatrix( './parameter_set/K=256_paraset/omega_frequencies_K=256.txt' );
				omega_GD = readmatrix( filepath_omega_GD );
				% eta = readmatrix( './parameter_set/K=256_paraset/weight_parameters_K=256.txt' );
				eta_GD = readmatrix( filepath_eta_GD );
				eta_GD_re = eta_GD( 1, : );
				eta_GD_im = eta_GD( 2, : );
				
				omega_AM = readmatrix( filepath_omega_AM );
				eta_AM = readmatrix( filepath_eta_AM );
				eta_AM_re = eta_AM( 1, : );
				eta_AM_im = eta_AM( 2, : );
				
				
				p_0_sample_store = randn( M_sp, dim ) * ( 1 / sqrt( beta ) );
				x_0_sample_store_AM = [];
				x_0_sample_store_GD = [];
				%%% Langevin dynamics sampling under the approximated potential with Adaptive Metropolis algorithm
				parfor j = 1 : Num_parallel_worker
					% seed_b = R * j^2;
					% rng( seed_b );
					x0_AM = randn( L_x_vectorize, dim );
					x0_GD = x0_AM;
					x0_sample_store_on_worker_AM = zeros( M_sp_on_worker, dim );
					x0_sample_store_on_worker_GD = zeros( M_sp_on_worker, dim );
					M_equil = 10^4;
					Delta_t_sample = 2;
					N_delta_t_sample = Delta_t_sample / h_x_sample;

					for n = 1 : 1 : M_equil
						Brownian_incre_equil_n = randn( L_x_vectorize, dim );
						delta_x_AM = -V_bar_re_prime( eta_AM_re, eta_AM_im, omega_AM, x0_AM, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) * h_x_equi + sqrt( h_x_equi ) * Brownian_incre_equil_n * sqrt( 2 / beta );
						x0_AM = x0_AM + delta_x_AM;
						delta_x_GD = -V_bar_re_prime( eta_AM_re, eta_AM_im, omega_AM, x0_GD, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) * h_x_equi + sqrt( h_x_equi ) * Brownian_incre_equil_n * sqrt( 2 / beta );
						x0_GD = x0_GD + delta_x_GD;
					end
					x_old_AM = x0_AM;
					x_old_GD = x0_GD;
					for n = 1 : 1 : M_sp_on_worker * N_delta_t_sample / L_x_vectorize
						Brownian_incre_sample_n = randn( L_x_vectorize, dim );
						delta_x_AM = -V_bar_re_prime( eta_AM_re, eta_AM_im, omega_AM, x_old_AM, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) * h_x_sample + sqrt( h_x_sample ) * Brownian_incre_sample_n * sqrt( 2 / beta );
						x_new_AM = x_old_AM + delta_x_AM;
						x_old_AM = x_new_AM;
						delta_x_GD = -V_bar_re_prime( eta_AM_re, eta_AM_im, omega_AM, x_old_GD, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) * h_x_sample + sqrt( h_x_sample ) * Brownian_incre_sample_n * sqrt( 2 / beta );
						x_new_GD = x_old_GD + delta_x_GD;
						x_old_GD = x_new_GD;
						
						if( mod( n, N_delta_t_sample ) == 0 )
							num_vectorize_x = n / N_delta_t_sample;
							x0_sample_store_on_worker_AM( ( num_vectorize_x - 1 ) * L_x_vectorize + 1 : num_vectorize_x * L_x_vectorize, : ) = x_new_AM;
							x0_sample_store_on_worker_GD( ( num_vectorize_x - 1 ) * L_x_vectorize + 1 : num_vectorize_x * L_x_vectorize, : ) = x_new_GD;
						end
					end
					% x_0_sample_store( ( ( j - 1 ) * M_sp_on_worker + 1 ) : j * M_sp_on_worker, : ) = x0_sample_store_on_worker;
					x_0_sample_store_AM = [ x_0_sample_store_AM; x0_sample_store_on_worker_AM ];
					x_0_sample_store_GD = [ x_0_sample_store_GD; x0_sample_store_on_worker_GD ];
					
				end
		
				% pr_z0_m = 1 / ( 2 * domain_bound )^2;
				pr_z0_m = 1;
				% Correlation function at time tau = 0
				Correlation_samples_AM_x( :, 1 ) = x_0_sample_store_AM( :, 1 ) .* x_0_sample_store_AM( :, 1 ) * pr_z0_m;
				Correlation_samples_AM_p( :, 1 ) = p_0_sample_store( :, 1 ) .* p_0_sample_store( :, 1 ) * pr_z0_m;
				Correlation_samples_GD_x( :, 1 ) = x_0_sample_store_GD( :, 1 ) .* x_0_sample_store_GD( :, 1 ) * pr_z0_m;
				Correlation_samples_GD_p( :, 1 ) = p_0_sample_store( :, 1 ) .* p_0_sample_store( :, 1 ) * pr_z0_m;
					
			
				% Correlation_samples_tau_posit_AM = zeros( M_sp, Num_tau_points - 1 );
				% Correlation_samples_tau_posit_GD = zeros( M_sp, Num_tau_points - 1 );
				Correlation_samples_tau_posit_AM_x = [];
				Correlation_samples_tau_posit_AM_p = [];
				Correlation_samples_tau_posit_GD_x = [];
				Correlation_samples_tau_posit_GD_p = [];
				% h_tau = 0.0025;
				h_tau = h_x_sample;
				n_steps_Delta_tau = round( Delta_tau_plot / h_tau );
				
			
				%%% Evaluating the approximation of the correlation function with Verlet method using the model trained by Adaptive Metropolis algorithm
				parfor( m = 1 : M_sp / L_x_vectorize )
					
					Correlation_sample_m_AM_x = zeros( L_x_vectorize, Num_tau_points - 1 );
					Correlation_sample_m_AM_p = zeros( L_x_vectorize, Num_tau_points - 1 );
					
					x_0_m = x_0_sample_store_AM( ( m - 1 ) * L_x_vectorize + 1 : m * L_x_vectorize, : );
					p_0_m = p_0_sample_store( ( m - 1 ) * L_x_vectorize + 1 : m * L_x_vectorize, : );

					x_old = x_0_m;
					p_old = p_0_m;
					x_n_new = x_0_m;
					p_n_new = p_0_m;
						
					for n = 1 : 1 : round( tau / h_tau )
						p_n_interme = p_old + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_AM_re, eta_AM_im, omega_AM, x_old, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) );
						x_n_new = x_old + h_tau * p_n_interme;
						p_n_new = p_n_interme + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_AM_re, eta_AM_im, omega_AM, x_n_new, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) );
						x_old = x_n_new;
						p_old = p_n_new;
						if( mod( n, n_steps_Delta_tau ) == 0 )
							num_plot_point = n / n_steps_Delta_tau;
							Correlation_sample_m_AM_x( :, num_plot_point ) = x_n_new( :, 1 ) .* x_0_m( :, 1 ) * pr_z0_m;
							Correlation_sample_m_AM_p( :, num_plot_point ) = p_n_new( :, 1 ) .* p_0_m( :, 1 ) * pr_z0_m;
						end
					end
					% Correlation_samples_tau_posit_AM( m, : ) = Correlation_sample_m_AM;
					Correlation_samples_tau_posit_AM_x = [ Correlation_samples_tau_posit_AM_x; Correlation_sample_m_AM_x ];
					Correlation_samples_tau_posit_AM_p = [ Correlation_samples_tau_posit_AM_p; Correlation_sample_m_AM_p ];
					
				end

				%%% Evaluating the approximation of the correlation function with Verlet method using the model trained by Gradient Descent algorithm
				parfor( m = 1 : M_sp / L_x_vectorize )
					
					Correlation_sample_m_GD_x = zeros( L_x_vectorize, Num_tau_points - 1 );
					Correlation_sample_m_GD_p = zeros( L_x_vectorize, Num_tau_points - 1 );
					
					x_0_m = x_0_sample_store_GD( ( m - 1 ) * L_x_vectorize + 1 : m * L_x_vectorize, : );
					p_0_m = p_0_sample_store( ( m - 1 ) * L_x_vectorize + 1 : m * L_x_vectorize, : );

					x_old = x_0_m;
					p_old = p_0_m;
					x_n_new = x_0_m;
					p_n_new = p_0_m;
						
					for n = 1 : 1 : round( tau / h_tau )
						p_n_interme = p_old + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_GD_re, eta_GD_im, omega_GD, x_old, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) );
						x_n_new = x_old + h_tau * p_n_interme;
						p_n_new = p_n_interme + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_GD_re, eta_GD_im, omega_GD, x_n_new, alpha, gamma, R_a, R_b, sigma_a, sigma_b ) );
						x_old = x_n_new;
						p_old = p_n_new;
						if( mod( n, n_steps_Delta_tau ) == 0 )
							num_plot_point = n / n_steps_Delta_tau;
							Correlation_sample_m_GD_x( :, num_plot_point ) = x_n_new( :, 1 ) .* x_0_m( :, 1 ) * pr_z0_m;
							Correlation_sample_m_GD_p( :, num_plot_point ) = p_n_new( :, 1 ) .* p_0_m( :, 1 ) * pr_z0_m;
						end
					end
					% Correlation_samples_tau_posit_AM( m, : ) = Correlation_sample_m_AM;
					Correlation_samples_tau_posit_GD_x = [ Correlation_samples_tau_posit_GD_x; Correlation_sample_m_GD_x ];
					Correlation_samples_tau_posit_GD_p = [ Correlation_samples_tau_posit_GD_p; Correlation_sample_m_GD_p ];
					
				end
				Correlation_samples_AM_x( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_AM_x;
				Correlation_samples_AM_p( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_AM_p;
				correlation_md_rep_R_AM_x = mean( Correlation_samples_AM_x, 1 );
				correlation_md_rep_R_AM_p = mean( Correlation_samples_AM_p, 1 );
				
				Correlation_samples_GD_x( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_GD_x;
				Correlation_samples_GD_p( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_GD_p;
				correlation_md_rep_R_GD_x = mean( Correlation_samples_GD_x, 1 );
				correlation_md_rep_R_GD_p = mean( Correlation_samples_GD_p, 1 );
				
				correlation_md_record_AM_x( R, : ) = correlation_md_rep_R_AM_x;
				correlation_md_record_AM_p( R, : ) = correlation_md_rep_R_AM_p;
				correlation_md_record_GD_x( R, : ) = correlation_md_rep_R_GD_x;
				correlation_md_record_GD_p( R, : ) = correlation_md_rep_R_GD_p;
				
			end
			
			t_MD = toc;
		    t_MD

			correlation_md_record_mean_AM_x = mean( correlation_md_record_AM_x, 1 );
			correlation_md_record_std_AM_x = std( correlation_md_record_AM_x, 0, 1 );
			correlation_md_record_mean_AM_p = mean( correlation_md_record_AM_p, 1 );
			correlation_md_record_std_AM_p = std( correlation_md_record_AM_p, 0, 1 );
			
			diff_MD_AM_exact_x = abs( correlation_curve_ref_x - correlation_md_record_mean_AM_x );
			L1_diff_rec_AM_x( i_a, 1 ) = mean( diff_MD_AM_exact_x );
			correlation_curve_ref_x_std = corr_curve_info_1( 2, 1 : Num_tau_points );
			L1_diff_rec_AM_x( i_a, 2 ) = mean( ( correlation_md_record_std_AM_x / sqrt( Num_replica ) + correlation_curve_ref_x_std / sqrt( 256 ) ) * 1.96 );
			
			diff_MD_AM_exact_p = abs( correlation_curve_ref_p - correlation_md_record_mean_AM_p );
			L1_diff_rec_AM_p( i_a, 1 ) = mean( diff_MD_AM_exact_p );
			correlation_curve_ref_p_std = corr_curve_info_1( 5, 1 : Num_tau_points );
			L1_diff_rec_AM_p( i_a, 2 ) = mean( ( correlation_md_record_std_AM_p / sqrt( Num_replica ) + correlation_curve_ref_p_std / sqrt( 256 ) ) * 1.96 );
			
			correlation_md_record_mean_GD_x = mean( correlation_md_record_GD_x, 1 );
			correlation_md_record_std_GD_x = std( correlation_md_record_GD_x, 0, 1 );
			correlation_md_record_mean_GD_p = mean( correlation_md_record_GD_p, 1 );
			correlation_md_record_std_GD_p = std( correlation_md_record_GD_p, 0, 1 );
			
			diff_MD_GD_exact_x = abs( correlation_curve_ref_x - correlation_md_record_mean_GD_x );
			L1_diff_rec_GD_x( i_a, 1 ) = mean( diff_MD_GD_exact_x );
			L1_diff_rec_GD_x( i_a, 2 ) = mean( ( correlation_md_record_std_GD_x / sqrt( Num_replica ) + correlation_curve_ref_x_std / sqrt( 256 ) ) * 1.96 );
			diff_MD_GD_exact_p = abs( correlation_curve_ref_p - correlation_md_record_mean_GD_p );
			L1_diff_rec_GD_p( i_a, 1 ) = mean( diff_MD_GD_exact_p );
			L1_diff_rec_GD_p( i_a, 2 ) = mean( ( correlation_md_record_std_GD_p / sqrt( Num_replica ) + correlation_curve_ref_p_std / sqrt( 256 ) ) * 1.96 );
			
			
			if( i_a == num_K_value )    % Plot the approximation of the correlation function for the largest K value
			
				% For the largest value of K, Plot the loaded model with the exact potential
				x1_points = ( -5 : 0.01 : 5 )';
				x2_points = ( -5 : 0.01 : 5 )';
				x_points = [ x1_points, x2_points ];
				num_x_points = length( x1_points );
				V_true_plot = zeros( num_x_points, 1 );
				V_int_bar_real_AM = zeros( num_x_points, 1 );
				V_int_bar_real_GD = zeros( num_x_points, 1 );
				% V_int_bar_real_prime = zeros( num_x_points, 1 );
				for i = 1 : 1 : num_x_points
					x_i = x_points( i, : );
					V_int_x_i_AM = 0;
					V_int_x_i_GD = 0;
					for k = 1 : 1 : length( eta_AM_re )
						omega_AM_k = omega_AM( :, k );
						inner_prod_AM = x_i * omega_AM_k;
						V_int_summand_AM = eta_AM_re( 1, k ) * cos( inner_prod_AM ) - eta_AM_im( 1, k ) * sin( inner_prod_AM );
						V_int_x_i_AM = V_int_x_i_AM + V_int_summand_AM;
						
						omega_GD_k = omega_GD( :, k );
						inner_prod_GD = x_i * omega_GD_k;
						V_int_summand_GD = eta_GD_re( 1, k ) * cos( inner_prod_GD ) - eta_GD_im( 1, k ) * sin( inner_prod_GD );
						V_int_x_i_GD = V_int_x_i_GD + V_int_summand_GD;
					end
					V_int_bar_real_AM( i, 1 ) = V_int_x_i_AM;
					V_int_bar_real_GD( i, 1 ) = V_int_x_i_GD;
					% V_true_plot( i, 1 ) = V_func( x_i, alpha, gamma, R_minus_2, R_minus_1 );
					V_true_plot( i, 1 ) = V_func( x_i, alpha, gamma, R_a, sigma_a );
				end
				figure( 1 )
				set( groot, 'defaultAxesTickLabelInterpreter', 'latex' ); 
				set( groot, 'defaultLegendInterpreter', 'latex' );
				plot( x1_points, V_true_plot );
				hold on
				plot( x1_points, V_int_bar_real_AM );
				plot( x1_points, V_int_bar_real_GD );
				hold off
				legend( 'Target potential $v(\mathbf{x})$', 'Adaptive Metropolis approximation $\bar{v}_{AM}(\mathbf{x})$', 'Gradient descent approximation $\bar{v}_{GD}(\mathbf{x})$' );
				
				
				fig2 = figure( 2 );
				tau_values = ( 0 : Delta_tau_plot : tau );
				set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
				set(groot, 'defaultLegendInterpreter','latex');
				
				hold on
				
				plot( tau_values, correlation_md_record_mean_AM_x );
				tau_values_2 = [ tau_values, fliplr( tau_values ) ];
				CI_upper_corr_curve = correlation_md_record_mean_AM_x + correlation_md_record_std_AM_x * 1.96 / sqrt( Num_replica );
				CI_lower_corr_curve = correlation_md_record_mean_AM_x - correlation_md_record_std_AM_x * 1.96 / sqrt( Num_replica );
				inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
				fill( tau_values_2, inBetween_statistic, 'cyan', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
				
				plot( tau_values, correlation_md_record_mean_GD_x );
				CI_upper_corr_curve = correlation_md_record_mean_GD_x + correlation_md_record_std_GD_x * 1.96 / sqrt( Num_replica );
				CI_lower_corr_curve = correlation_md_record_mean_GD_x - correlation_md_record_std_GD_x * 1.96 / sqrt( Num_replica );
				inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
				fill( tau_values_2, inBetween_statistic, 'magenta', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
				
				plot( tau_values, correlation_curve_ref_x );
				CI_upper_corr_curve = correlation_curve_ref_x + correlation_curve_ref_x_std * 1.96 / sqrt( Num_replica );
				CI_lower_corr_curve = correlation_curve_ref_x - correlation_curve_ref_x_std * 1.96 / sqrt( Num_replica );
				inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
				fill( tau_values_2, inBetween_statistic, 'green', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
				
				hold off
				legend( 'MD with Adaptive Metropolis trained $\bar{v}_{AM}(x)$', 'CI of Adaptive Metropolis MD correlation', 'MD with Gradient descent trained $\bar{v}_{GD}(x)$', 'CI of Gradient Descent MD correlation', 'Reference curve', 'CI of the reference correlation curve' );
				xlabel( 'correlation time $\tau$' )
				ylabel( '$\langle x_1(\tau),x_1(0)\rangle$-correlation' )
				
				fig3 = figure( 3 );
				tau_values = ( 0 : Delta_tau_plot : tau );
				set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
				set(groot, 'defaultLegendInterpreter','latex');

				hold on
				
				plot( tau_values, correlation_md_record_mean_AM_p );
				tau_values_2 = [ tau_values, fliplr( tau_values ) ];
				CI_upper_corr_curve = correlation_md_record_mean_AM_p + correlation_md_record_std_AM_p * 1.96 / sqrt( Num_replica );
				CI_lower_corr_curve = correlation_md_record_mean_AM_p - correlation_md_record_std_AM_p * 1.96 / sqrt( Num_replica );
				inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
				fill( tau_values_2, inBetween_statistic, 'cyan', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
				
				plot( tau_values, correlation_md_record_mean_GD_p );
				CI_upper_corr_curve = correlation_md_record_mean_GD_p + correlation_md_record_std_GD_p * 1.96 / sqrt( Num_replica );
				CI_lower_corr_curve = correlation_md_record_mean_GD_p - correlation_md_record_std_GD_p * 1.96 / sqrt( Num_replica );
				inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
				fill( tau_values_2, inBetween_statistic, 'magenta', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
				
				plot( tau_values, correlation_curve_ref_p );
				CI_upper_corr_curve = correlation_curve_ref_p + correlation_curve_ref_p_std * 1.96 / sqrt( Num_replica );
				CI_lower_corr_curve = correlation_curve_ref_p - correlation_curve_ref_p_std * 1.96 / sqrt( Num_replica );
				inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
				fill( tau_values_2, inBetween_statistic, 'green', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
				
				hold off
				legend( 'MD with Adaptive Metropolis trained $\bar{v}_{AM}(x)$', 'CI of Adaptive Metropolis MD correlation', 'MD with Gradient descent trained $\bar{v}_{GD}(x)$', 'CI of Gradient Descent MD correlation', 'Reference curve', 'CI of the reference correlation curve' );
				xlabel( 'correlation time $\tau$' )
				ylabel( '$\langle p_1(\tau),p_1(0)\rangle$-correlation' )
			end
		
		end
		poolobj = gcp( 'nocreate' );
		delete( poolobj )
		
		
		% profile viewer
		%{
		correlation_exact = zeros( Num_tau_points, 1 );
		for i = 1 : 1 : Num_tau_points
			tau_i = tau_values( i, 1 );
			correlation_exact( i, 1 ) = cos( tau_i ) / beta;
		end
		%}
		
		
		fig4 = figure( 4 );
		% plot( K_values, L1_diff_rec_AM_x( :, 1 ), '-*' );
		hold on
		errorbar( K_values, L1_diff_rec_AM_x( :, 1 ), L1_diff_rec_AM_x( :, 2 ), '-*' );
		% plot( K_values, L1_diff_rec_GD_x( :, 1 ), '-o' );
		errorbar( K_values, L1_diff_rec_GD_x( :, 1 ), L1_diff_rec_GD_x( :, 2 ), '-o' );
		% loglog( K_values_sgd, Errors_K_rec_3, '-+' )
		plot( K_values, 0.5 * K_values.^(-1/2) );
		hold off
		% legend( 'loss for adaptive Metropolis algorithm', 'loss for Gradient Descent on frequency', 'loss for SGD on frequency and amplitude', 'reference $\mathcal{O}(K^{-1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
		legend( '$L_1$ difference of adaptive Metropolis algorithm', '$L_1$ difference of  Gradient Descent on frequency', 'reference $\mathcal{O}(K^{-1/2})$', 'fontsize', 20, 'interpreter', 'latex' ); 
		% legend( 'loss for adaptive Metropolis algorithm', 'reference $\mathcal{O}(K^{-1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
		xlabel( '$K$', 'fontsize', 20, 'interpreter', 'latex' );
		ylabel( '$L_1$ difference of correlation function', 'fontsize', 20, 'interpreter', 'latex' );
		% title( 'Log-log plot of testing loss, target function $f(x)=e^{-\frac{|x|^{2}}{2}}\,\mathrm{Si}(\frac{x}{a})$, $a=1$, $J=10^3$, $d=1$', 'fontsize', 20, 'interpreter', 'latex' )
		% title( 'Log-log plot of testing loss, target function $v_1(x)=\frac{|x_1|^2}{2}\,\chi_1(|x|)$, $J=10^3$, $d=2$', 'fontsize', 20, 'interpreter', 'latex' )
		title_string = sprintf('Log-log plot of $L_1$ difference in the $\langle x_1(0), x_1(\tau) \rangle$ auto-correlation function, $J=10^{5}$, $d=%d$', dim);
		title( title_string, 'fontsize', 20, 'interpreter', 'latex' )
		set( gca, 'YScale', 'log');
		set( gca, 'XScale', 'log');
		
		fig5 = figure( 5 );
		hold on
		errorbar( K_values, L1_diff_rec_AM_p( :, 1 ), L1_diff_rec_AM_p( :, 2 ), '-*' );
		errorbar( K_values, L1_diff_rec_GD_p( :, 1 ), L1_diff_rec_GD_p( :, 2 ), '-o' );
		% loglog( K_values_sgd, Errors_K_rec_3, '-+' )
		plot( K_values, 0.5 * K_values.^(-1/2) );
		hold off
		legend( '$L_1$ difference of adaptive Metropolis algorithm', '$L_1$ difference of  Gradient Descent on frequency', 'reference $\mathcal{O}(K^{-1/2})$', 'fontsize', 20, 'interpreter', 'latex' ); 
		xlabel( '$K$', 'fontsize', 20, 'interpreter', 'latex' );
		ylabel( '$L_1$ difference of correlation function', 'fontsize', 20, 'interpreter', 'latex' );
		title_string = sprintf('Log-log plot of $L_1$ difference in the $\langle p_1(0), p_1(\tau) \rangle$ auto-correlation function, $J=10^{5}$, $d=%d$', dim);
		title( title_string, 'fontsize', 20, 'interpreter', 'latex' )
		set( gca, 'YScale', 'log');
		set( gca, 'XScale', 'log');


        folderName = 'Figures_save';
		if ~exist( folderName, 'dir' )
			mkdir( folderName );
		end
		corr_figure_save_filename_2 = sprintf( 'correlation_md_x_AM_h=%6.5f_M_sp=%d_N_rep=%d.fig', h_tau, M_sp, Num_replica );
		saveas( fig2, fullfile( folderName, corr_figure_save_filename_2 ) ); % Using figure number
		
		corr_figure_save_filename_3 = sprintf( 'correlation_md_p_AM_h=%6.5f_M_sp=%d_N_rep=%d.fig', h_tau, M_sp, Num_replica );
		saveas( fig3, fullfile( folderName, corr_figure_save_filename_3 ) ); % Using figure number

        corr_figure_save_filename_4 = sprintf( 'correlation_L1_diff_x_h=%6.5f_M_sp=%d_N_rep=%d.fig', h_tau, M_sp, Num_replica );
		saveas( fig4, fullfile( folderName, corr_figure_save_filename_4 ) ); % Using figure number
        
        corr_figure_save_filename_5 = sprintf( 'correlation_L1_diff_p_h=%6.5f_M_sp=%d_N_rep=%d.fig', h_tau, M_sp, Num_replica );
		saveas( fig5, fullfile( folderName, corr_figure_save_filename_5 ) ); % Using figure number
		
		corr_curve_info = 1;

	end
	
	
	
	
	
	
	
	
	
	
	
	