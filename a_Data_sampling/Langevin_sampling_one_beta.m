    
	% clear all
    % close all
	
	function Langevin_sampling_test = Langevin_sampling_one_beta( beta, dim, alpha, gamma, R_a, R_c, sigma_a, sigma_c, h_x_sample, M_sp, Num_replica, Num_parallel_worker, L_x_vectorize )
	
		% dim = 2;
		% beta = 0.2;
		% alpha = sqrt( 2 );
		% gamma = 2;
		
		% M_sp = 10^6;
		h_x_equi = 0.1;
		% h_x_sample = 0.01;
		x_0_sample_store = [];
		
		parpool( Num_parallel_worker );
		M_sp_on_worker = M_sp / Num_parallel_worker;
		
		tic
		
		for Rep = 1 : 1 : Num_replica
            % seed_a = Rep^2;
            % rng( seed_a );
			x_0_sample_store_Rep = [];
			parfor j = 1 : Num_parallel_worker
				% seed_b = R * j^2;
				% rng( seed_b );
				x0 = randn( L_x_vectorize, dim );
				x0_sample_store_on_worker = zeros( M_sp_on_worker, dim );
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
				x_0_sample_store_Rep = [ x_0_sample_store_Rep; x0_sample_store_on_worker ];
			end
			x_0_sample_store = [ x_0_sample_store; x_0_sample_store_Rep ];
		end
		
		t_sampling = toc;
		t_sampling

		
		V_x0_sample_store = zeros( M_sp * Num_replica, 1 );
		V_prime_x0_sample_store = zeros( M_sp * Num_replica, dim );
		parfor m = 1 : 1 : M_sp * Num_replica
			x_m = x_0_sample_store( m, : );
			% V_x0_sample_store( m, 1 ) = V_func( x_m, alpha, gamma, R_minus_2, R_minus_1 );
			% V_prime_x0_sample_store( m, : ) = V_prime( x_m, alpha, gamma, R_minus_2, R_minus_1 );
			V_x0_sample_store( m, 1 ) = V_1_func( x_m, alpha, gamma, R_a, R_c, sigma_a, sigma_c );
			V_prime_x0_sample_store( m, : ) = V_1_prime( x_m, alpha, gamma, R_a, R_c, sigma_a, sigma_c );
		end
		
		poolobj = gcp( 'nocreate' );
		delete( poolobj );
		
		figure( 3 )
		plot( x_0_sample_store( :, 1 ), x_0_sample_store( :, 2 ), '.' );
		
		
		if ( ~exist( 'Sample_store', 'dir' ) )
			mkdir( 'Sample_store' );
        end
        Delta_t_sample = 2;
		sample_matrix = [ x_0_sample_store, V_x0_sample_store, V_prime_x0_sample_store ];
		filename = "Sample_store\Langevin_sample_store_V_general_one_beta" + "_dim=" + dim + "_Msp=" + M_sp + "_Num_rep=" + Num_replica + "_Delta_t_sample=" + Delta_t_sample + "_dt=" + h_x_sample + ".csv";
		writematrix( sample_matrix, filename )


		figure( 4 )
		histogram( x_0_sample_store( :, 1 ) );
	
		Langevin_sampling_test = 1;
	end