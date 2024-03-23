    
	function corr_curve_info = corr_function_curve( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, Num_parallel_worker, L_x_vectorize )
		
		beta = 1;
		dim = 2;
		x_corr = 1;
		alpha = sqrt( 2 );
		gamma = 2;
		
		R_a = 2;    % R_b = 3;
		R_c = 4;    
		sigma_a = 1;
		
		
		% tau = 2;
		% Delta_tau_plot = 0.1;
		Num_tau_points = round( tau / Delta_tau_plot ) + 1;
		% M_sp = 10^5;    % For 10 CPU cores
		% M_sp = 2^12;    % For 8 CPU cores
		Correlation_samples_x = zeros( M_sp, Num_tau_points );
		Correlation_samples_p = zeros( M_sp, Num_tau_points );
		Correlation_samples_x_2h = zeros( M_sp, Num_tau_points );
		Correlation_samples_p_2h = zeros( M_sp, Num_tau_points );
		
		h_x_equi = 0.1;
		% h_x_sample = 0.025;
		
		% Num_parallel_worker = 8;
		parpool( Num_parallel_worker );
		M_sp_on_worker = M_sp / Num_parallel_worker;
		% L_x_vectorize = 1;
		
		tic
		
		% Num_replica = 32;
		correlation_md_record_x = zeros( Num_replica, Num_tau_points );
		correlation_md_record_p = zeros( Num_replica, Num_tau_points );
		correlation_md_record_x_2h = zeros( Num_replica, Num_tau_points );
		correlation_md_record_p_2h = zeros( Num_replica, Num_tau_points );
		
		% seed = 428;
		% rng( seed );
		for R = 1 : 1 : Num_replica
            % seed_a = R^2;
            % rng( seed_a );
			p_0_sample_store = randn( M_sp, dim ) * ( 1 / sqrt( beta ) );
			x_0_sample_store = [];
			parfor j = 1 : Num_parallel_worker
                % seed_b = R * j^2;
                % rng( seed_b );
				x0 = randn( L_x_vectorize, dim );
				x0_sample_store_on_worker = zeros( M_sp_on_worker, dim );
				% x0_sample_store_on_worker_V_true = zeros( M_sp_on_worker, dim );
				M_equil = 10^4;
				% Brownian_incre_equil = randn( L_x_vectorize, dim, M_equil );
				Delta_t_sample = 2;
				N_delta_t_sample = Delta_t_sample / h_x_sample;
				% Brownian_incre_sample = randn( N_delta_t_sample * M_sp_on_worker, dim );
				for n = 1 : 1 : M_equil
					Brownian_incre_equil_n = randn( L_x_vectorize, dim );
					delta_x = -V_prime( x0, alpha, gamma, R_a, sigma_a ) * h_x_equi + sqrt( h_x_equi ) * Brownian_incre_equil_n * sqrt( 2 / beta );
					x0 = x0 + delta_x;
				end
				x_old = x0;
				for n = 1 : 1 : M_sp_on_worker * N_delta_t_sample / L_x_vectorize
					Brownian_incre_sample_n = randn( L_x_vectorize, dim );
					% delta_x = -V_prime( x_old, alpha, gamma, R_minus_2, alpha_chi ) * h_x_sample + sqrt( h_x_sample ) * Brownian_incre_sample( n, : ) * sqrt( 2 / beta );
					delta_x = -V_prime( x_old, alpha, gamma, R_a, sigma_a ) * h_x_sample + sqrt( h_x_sample ) * Brownian_incre_sample_n * sqrt( 2 / beta );
					x_new = x_old + delta_x;
					x_old = x_new;
					if( mod( n, N_delta_t_sample ) == 0 )
						num_vectorize_x = n / N_delta_t_sample;
						% x0_sample_store_on_worker( n / N_delta_t_sample, : ) = x_new;
						x0_sample_store_on_worker( ( num_vectorize_x - 1 ) * L_x_vectorize + 1 : num_vectorize_x * L_x_vectorize, : ) = x_new;
					end
				end
				
				% x_0_sample_store( ( ( j - 1 ) * M_sp_on_worker + 1 ) : j * M_sp_on_worker, : ) = x0_sample_store_on_worker;
				x_0_sample_store = [ x_0_sample_store; x0_sample_store_on_worker ];
			end
			
			% pr_z0_m = 1 / ( 2 * domain_bound )^2;
			pr_z0_m = 1;
			% Correlation function at time tau = 0
			Correlation_samples_x( :, 1 ) = x_0_sample_store( :, 1 ) .* x_0_sample_store( :, 1 ) * pr_z0_m;
			Correlation_samples_p( :, 1 ) = p_0_sample_store( :, 1 ) .* p_0_sample_store( :, 1 ) * pr_z0_m;
			
			Correlation_samples_x_2h( :, 1 ) = x_0_sample_store( :, 1 ) .* x_0_sample_store( :, 1 ) * pr_z0_m;
			Correlation_samples_p_2h( :, 1 ) = p_0_sample_store( :, 1 ) .* p_0_sample_store( :, 1 ) * pr_z0_m;

			% t_sampling = toc
					
			% Correlation_samples_tau_posit_x = zeros( M_sp, Num_tau_points - 1 );
            % Correlation_samples_tau_posit_p = zeros( M_sp, Num_tau_points - 1 );
			% Correlation_samples_tau_posit_x_2h = zeros( M_sp, Num_tau_points - 1 );
            % Correlation_samples_tau_posit_p_2h = zeros( M_sp, Num_tau_points - 1 );
			Correlation_samples_tau_posit_x = [];
            Correlation_samples_tau_posit_p = [];
			Correlation_samples_tau_posit_x_2h = [];
            Correlation_samples_tau_posit_p_2h = [];
			% h_tau = 0.005;
			h_tau = h_x_sample;
			h_tau_double = 2 * h_tau;
			n_steps_Delta_tau = round( Delta_tau_plot / h_tau );
			n_steps_Delta_tau_2h = round( Delta_tau_plot / h_tau_double );
			
			% tic 
			parfor( m = 1 : M_sp / L_x_vectorize )
				
				Correlation_sample_m_x = zeros( L_x_vectorize, Num_tau_points - 1 );
                Correlation_sample_m_p = zeros( L_x_vectorize, Num_tau_points - 1 );
				Correlation_sample_m_x_2h = zeros( L_x_vectorize, Num_tau_points - 1 );
                Correlation_sample_m_p_2h = zeros( L_x_vectorize, Num_tau_points - 1 );
				
				% x_0_m = x_0_sample_store( m, : );
				% p_0_m = p_0_sample_store( m, : );
				x_0_m = x_0_sample_store( ( m - 1 ) * L_x_vectorize + 1 : m * L_x_vectorize, : );
				p_0_m = p_0_sample_store( ( m - 1 ) * L_x_vectorize + 1 : m * L_x_vectorize, : );

				x_old = x_0_m;
				p_old = p_0_m;
				x_n_new = x_0_m;
				p_n_new = p_0_m;
					
				for n = 1 : 1 : round( tau / h_tau )
					% p_n_interme = p_old + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_re, eta_im, omega, x0, alpha, gamma, R_minus_2, R_minus_1 ) );
					p_n_interme = p_old + ( h_tau / 2 ) * ( -V_prime( x_old, alpha, gamma, R_a, sigma_a ) );
					x_n_new = x_old + h_tau * p_n_interme;
					p_n_new = p_n_interme + ( h_tau / 2 ) * ( -V_prime( x_n_new, alpha, gamma, R_a, sigma_a ) );
					% p_n_new = p_n_interme + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_re, eta_im, omega, x0, alpha, gamma, R_minus_2, R_minus_1 ) );
					x_old = x_n_new;
					p_old = p_n_new;
					if( mod( n, n_steps_Delta_tau ) == 0 )
						num_plot_point = n / n_steps_Delta_tau;
						% if( x_corr == 1 )
						Correlation_sample_m_x( :, num_plot_point ) = x_n_new( :, 1 ) .* x_0_m( :, 1 ) * pr_z0_m;
						% elseif( x_corr == 0 )
						Correlation_sample_m_p( :, num_plot_point ) = p_n_new( :, 1 ) .* p_0_m( :, 1 ) * pr_z0_m;
						% end
					end
				end
				% Correlation_samples_tau_posit_x( m, : ) = Correlation_sample_m_x;
				% Correlation_samples_tau_posit_p( m, : ) = Correlation_sample_m_p;
				Correlation_samples_tau_posit_x = [ Correlation_samples_tau_posit_x; Correlation_sample_m_x ];
				Correlation_samples_tau_posit_p = [ Correlation_samples_tau_posit_p; Correlation_sample_m_p ];
				
				x_old = x_0_m;
				p_old = p_0_m;
				x_n_new = x_0_m;
				p_n_new = p_0_m;
				for n = 1 : 1 : round( tau / h_tau_double )
					% p_n_interme = p_old + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_re, eta_im, omega, x0, alpha, gamma, R_minus_2, R_minus_1 ) );
					p_n_interme = p_old + ( h_tau_double / 2 ) * ( -V_prime( x_old, alpha, gamma, R_a, sigma_a ) );
					x_n_new = x_old + h_tau_double * p_n_interme;
					p_n_new = p_n_interme + ( h_tau_double / 2 ) * ( -V_prime( x_n_new, alpha, gamma, R_a, sigma_a ) );
					% p_n_new = p_n_interme + ( h_tau / 2 ) * ( -V_bar_re_prime( eta_re, eta_im, omega, x0, alpha, gamma, R_minus_2, R_minus_1 ) );
					x_old = x_n_new;
					p_old = p_n_new;
					if( mod( n, n_steps_Delta_tau_2h ) == 0 )
						num_plot_point = n / n_steps_Delta_tau_2h;
						% if( x_corr == 1 )
						Correlation_sample_m_x_2h( :, num_plot_point ) = x_n_new( :, 1 ) .* x_0_m( :, 1 ) * pr_z0_m;
						% elseif( x_corr == 0 )
						Correlation_sample_m_p_2h( :, num_plot_point ) = p_n_new( :, 1 ) .* p_0_m( :, 1 ) * pr_z0_m;
						% end
					end
				end

				% Correlation_samples_tau_posit_x_2h( m, : ) = Correlation_sample_m_x_2h;
				% Correlation_samples_tau_posit_p_2h( m, : ) = Correlation_sample_m_p_2h;
				Correlation_samples_tau_posit_x_2h = [ Correlation_samples_tau_posit_x_2h; Correlation_sample_m_x_2h ];
				Correlation_samples_tau_posit_p_2h = [ Correlation_samples_tau_posit_p_2h; Correlation_sample_m_p_2h ];
				
            end
            Correlation_samples_x( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_x;
			Correlation_samples_p( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_p;
			
			Correlation_samples_x_2h( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_x_2h;
			Correlation_samples_p_2h( :, 2 : Num_tau_points ) = Correlation_samples_tau_posit_p_2h;
		
			correlation_md_rep_R_x = mean( Correlation_samples_x, 1 );
			correlation_md_rep_R_p = mean( Correlation_samples_p, 1 );
			correlation_md_rep_R_x_2h = mean( Correlation_samples_x_2h, 1 );
			correlation_md_rep_R_p_2h = mean( Correlation_samples_p_2h, 1 );
			
			correlation_md_record_x( R, : ) = correlation_md_rep_R_x;
			correlation_md_record_p( R, : ) = correlation_md_rep_R_p;
			correlation_md_record_x_2h( R, : ) = correlation_md_rep_R_x_2h;
			correlation_md_record_p_2h( R, : ) = correlation_md_rep_R_p_2h;
		end
		
		
		poolobj = gcp( 'nocreate' );
		delete( poolobj )
		
		t_MD = toc;
		t_MD
		
		correlation_md_record_mean_x = mean( correlation_md_record_x, 1 );
		correlation_md_record_std_x = std( correlation_md_record_x, 0, 1 );
		correlation_md_record_mean_p = mean( correlation_md_record_p, 1 );
		correlation_md_record_std_p = std( correlation_md_record_p, 0, 1 );
		
		correlation_md_record_mean_x_2h = mean( correlation_md_record_x_2h, 1 );
		correlation_md_record_mean_p_2h = mean( correlation_md_record_p_2h, 1 );
		
		corr_curve_info = [ correlation_md_record_mean_x; correlation_md_record_std_x; correlation_md_record_mean_x_2h;
							correlation_md_record_mean_p; correlation_md_record_std_p; correlation_md_record_mean_p_2h ];
		
	
	end
		
		
		
		
		
		
		
	
	
	
	