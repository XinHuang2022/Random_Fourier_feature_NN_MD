    
	% function corr_curve_info = test_approx_V_md_corr( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, Num_parallel_worker )
		
		Num_replica = 32;
		% C_const_bound = 200;
        C_const_bound = 0.1;
		
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
		% K_values = [ 8, 16, 32, 64, 128, 256, 512 ];
		K_values = [ 8, 16, 32, 64, 128, 256 ];
		% K_values = 64;
		num_K_value = length( K_values );
		
		alpha_3_record_AM = zeros( num_K_value, 2 );    % first row for storing the mean of the L1-diff with R replica, second row for the std of the L1-diff with R replica.
		alpha_3_record_GD = zeros( num_K_value, 2 );
		
        tic
		for i_a = 1 : 1 : num_K_value
		
			K = K_values( 1, i_a );
			
			alpha_3_samples_replica_AM = zeros( Num_replica, 1 );
			alpha_3_samples_replica_GD = zeros( Num_replica, 1 );
			
			
			for R = 1 : Num_replica
				
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
				
				sum_indicator_norm_AM = 0;
				sum_penalty_term_AM = 0;
				sum_indicator_norm_GD = 0;
				sum_penalty_term_GD = 0;
				
				for k = 1 : 1 : K
					eta_k_AM_real = eta_AM_re( 1, k );
					eta_k_AM_imag = eta_AM_im( 1, k );
					omega_k_AM = omega_AM( :, k );
					
					summand_indicator_norm = sqrt( eta_k_AM_real^2 + eta_k_AM_imag^2 ) * ( 1 + norm( omega_k_AM ) + norm( omega_k_AM * omega_k_AM', "fro" ) );
					sum_indicator_norm_AM = sum_indicator_norm_AM + summand_indicator_norm;
					
					summand_penalty_term = sqrt( eta_k_AM_real^2 + eta_k_AM_imag^2 ) * ( 1 + norm( omega_k_AM ) );
					% sum_penalty_term_AM = sum_penalty_term_AM + summand_penalty_term;
                    sum_penalty_term_AM = sum_indicator_norm_AM;
					
					eta_k_GD_real = eta_GD_re( 1, k );
					eta_k_GD_imag = eta_GD_im( 1, k );
					omega_k_GD = omega_GD( :, k );
					
					summand_indicator_norm = sqrt( eta_k_GD_real^2 + eta_k_GD_imag^2 ) * ( 1 + norm( omega_k_GD ) + norm( omega_k_GD * omega_k_GD', "fro" ) );
                    sum_indicator_norm_GD = sum_indicator_norm_GD + summand_indicator_norm;
					
					summand_penalty_term = sqrt( eta_k_GD_real^2 + eta_k_GD_imag^2 ) * ( 1 + norm( omega_k_GD ) );
					% sum_penalty_term_GD = sum_penalty_term_GD + summand_penalty_term;
                    sum_penalty_term_GD = sum_indicator_norm_GD;
				end
				
				if( sum_indicator_norm_AM > C_const_bound )
					% alpha_3_samples_replica_AM( R, 1 ) = sum_penalty_term_AM^2 + 1;
                    alpha_3_samples_replica_AM( R, 1 ) = sum_penalty_term_AM;
				end
				if( sum_indicator_norm_GD > C_const_bound )
					% alpha_3_samples_replica_GD( R, 1 ) = sum_penalty_term_GD^2 + 1;
                    alpha_3_samples_replica_GD( R, 1 ) = sum_penalty_term_GD;
				end
				
			end
			
			alpha_3_record_AM( i_a, 1 ) = mean( alpha_3_samples_replica_AM );
			alpha_3_record_AM( i_a, 2 ) = std( alpha_3_samples_replica_AM );
			alpha_3_record_GD( i_a, 1 ) = mean( alpha_3_samples_replica_GD );
			alpha_3_record_GD( i_a, 2 ) = std( alpha_3_samples_replica_GD );
			
			

		end
		t_MD = toc;
		t_MD
		
		fig9 = figure( 9 );
		hold on
		errorbar( K_values, alpha_3_record_AM( :, 1 ), alpha_3_record_AM( :, 2 ), '-*' );
		errorbar( K_values, alpha_3_record_GD( :, 1 ), alpha_3_record_GD( :, 2 ), '-o' );
		plot( K_values, 0.5 * K_values.^(1) );
		hold off
		legend( '$\alpha_3$ penalty of adaptive Metropolis algorithm', '$\alpha_3$ penalty of  Gradient Descent on frequency', 'reference $\mathcal{O}(K^{1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
		xlabel( '$K$', 'fontsize', 20, 'interpreter', 'latex' );
		ylabel( '$\alpha_3$ penalty term', 'fontsize', 20, 'interpreter', 'latex' );
		title_string = sprintf('Log-log plot of $\alpha_3$ penalty term, $J=10^{4}$, $d=%d$', dim);
		title( title_string, 'fontsize', 20, 'interpreter', 'latex' )
		set( gca, 'YScale', 'log');
		set( gca, 'XScale', 'log');
		
		
	
	
	
	
	
	
	
	
	
	
	