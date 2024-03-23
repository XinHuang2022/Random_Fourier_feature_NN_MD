    

	function Langevin_sampling_test = Langevin_sampling_mix_beta( beta_1, beta_2, dim, alpha, gamma, R_a, R_c, sigma_a, sigma_c, h_x_sample, M_sp, Num_replica, Num_parallel_worker, N, Training_data_store_path )
	
		h_x_equi = 0.1;					% h_x_equi is the time step length before equilibrium in the Langevin dynamics sampling
		L_x_vectorize = 16;				% L_x_vectorize is the length for vectorization to speed up the generating of samples
		Delta_t_sample = 2;             % Delta_t_sample is the time difference between two samples with overdamped Langevin dynamics
        
        x_0_sample_store = [];
		
		parpool( Num_parallel_worker );
		M_sp_on_worker = M_sp / Num_parallel_worker;
		
		tic
		
		for Rep = 1 : 1 : ( Num_replica / 2 )
            % seed_a = Rep^2;
            % rng( seed_a );
			x_0_sample_store_Rep = [];
			parfor j = 1 : Num_parallel_worker
				% seed_b = R * j^2;
				% rng( seed_b );
				x0 = randn( L_x_vectorize, dim );
				x0_sample_store_on_worker = zeros( M_sp_on_worker, dim );
				M_equil = 10^4;
				
				N_delta_t_sample = Delta_t_sample / h_x_sample;
				for n = 1 : 1 : M_equil
					Brownian_incre_equil_n = randn( L_x_vectorize, dim );
					delta_x = -V_prime( x0, alpha, gamma, R_a, sigma_a ) * h_x_equi + sqrt( h_x_equi ) * Brownian_incre_equil_n * sqrt( 2 / beta_1 );
					x0 = x0 + delta_x;
				end
				x_old = x0;
				for n = 1 : 1 : M_sp_on_worker * N_delta_t_sample / L_x_vectorize
					Brownian_incre_sample_n = randn( L_x_vectorize, dim );
					delta_x = -V_prime( x_old, alpha, gamma, R_a, sigma_a ) * h_x_sample + sqrt( h_x_sample ) * Brownian_incre_sample_n * sqrt( 2 / beta_1 );
					x_new = x_old + delta_x;
					x_old = x_new;
					if( mod( n, N_delta_t_sample ) == 0 )
						num_vectorize_x = n / N_delta_t_sample;
						x0_sample_store_on_worker( ( num_vectorize_x - 1 ) * L_x_vectorize + 1 : num_vectorize_x * L_x_vectorize, : ) = x_new;
					end
				end
				x_0_sample_store_Rep = [ x_0_sample_store_Rep; x0_sample_store_on_worker ];
			end
			x_0_sample_store = [ x_0_sample_store; x_0_sample_store_Rep ];
		end
		for Rep = 1 : 1 : ( Num_replica / 2 )
            % seed_a = Rep^2;
            % rng( seed_a );
			x_0_sample_store_Rep = [];
			parfor j = 1 : Num_parallel_worker
				% seed_b = R * j^2;
				% rng( seed_b );
				x0 = randn( L_x_vectorize, dim );
				x0_sample_store_on_worker = zeros( M_sp_on_worker, dim );
				M_equil = 10^4;
				N_delta_t_sample = Delta_t_sample / h_x_sample;
				for n = 1 : 1 : M_equil
					Brownian_incre_equil_n = randn( L_x_vectorize, dim );
					delta_x = -V_prime( x0, alpha, gamma, R_a, sigma_a ) * h_x_equi + sqrt( h_x_equi ) * Brownian_incre_equil_n * sqrt( 2 / beta_2 );
					x0 = x0 + delta_x;
				end
				x_old = x0;
				for n = 1 : 1 : M_sp_on_worker * N_delta_t_sample / L_x_vectorize
					Brownian_incre_sample_n = randn( L_x_vectorize, dim );
					delta_x = -V_prime( x_old, alpha, gamma, R_a, sigma_a ) * h_x_sample + sqrt( h_x_sample ) * Brownian_incre_sample_n * sqrt( 2 / beta_2 );
					x_new = x_old + delta_x;
					x_old = x_new;
					if( mod( n, N_delta_t_sample ) == 0 )
						num_vectorize_x = n / N_delta_t_sample;
						x0_sample_store_on_worker( ( num_vectorize_x - 1 ) * L_x_vectorize + 1 : num_vectorize_x * L_x_vectorize, : ) = x_new;
					end
				end
				x_0_sample_store_Rep = [ x_0_sample_store_Rep; x0_sample_store_on_worker ];
			end
			x_0_sample_store = [ x_0_sample_store; x_0_sample_store_Rep ];
		end
		
		t_sampling = toc;
		t_sampling

		
		V_x0_sample_store = zeros( M_sp * Num_replica, 1 );
		V_prime_x0_sample_store = zeros( M_sp * Num_replica, dim );
		parfor m = 1 : M_sp * Num_replica
			x_m = x_0_sample_store( m, : );
			V_x0_sample_store( m, 1 ) = V_1_func( x_m, alpha, gamma, R_a, R_c, sigma_a, sigma_c );
			V_prime_x0_sample_store( m, : ) = V_1_prime( x_m, alpha, gamma, R_a, R_c, sigma_a, sigma_c );
		end
		
		poolobj = gcp( 'nocreate' );
		delete( poolobj );
		
		Plot_dots_number = 10^6;
		if( M_sp < Plot_dots_number )
			fig1 = figure( 1 );
			plot( x_0_sample_store( :, 1 ), x_0_sample_store( :, 2 ), '.' );
			Fig_FolderName = fullfile( 'a_Data_sampling_pipeline', 'Figures_save_data_sampling' );
			if ~exist( Fig_FolderName, 'dir' )			
				mkdir( Fig_FolderName );
			end
			fig1_name = "Sampled_data_plot_mix_beta" + "_dim=" + dim + "_Msp=" + M_sp + "_Num_rep=" + Num_replica + "_Delta_t_sample=" + Delta_t_sample + "_dt=" + h_x_sample + ".fig";
			saveas( fig1, fullfile( Fig_FolderName, fig1_name ) ); 
			
			fig2 = figure( 2 );
			histogram( x_0_sample_store( :, 1 ) );
			fig2_name = "Sampled_data_histogram_mix_beta" + "_dim=" + dim + "_Msp=" + M_sp + "_Num_rep=" + Num_replica + "_Delta_t_sample=" + Delta_t_sample + "_dt=" + h_x_sample + ".fig";
			saveas( fig2, fullfile( Fig_FolderName, fig2_name ) );
		end
		
		sample_matrix = [ x_0_sample_store, V_x0_sample_store, V_prime_x0_sample_store ];
		% Reshuffle on the data set with mixed beta values
		[ N_row, N_col ] = size( sample_matrix );
		rand_indices = randperm( N_row );
		data_set = zeros( N_row, N_col );
		data_set( rand_indices, : )  = sample_matrix;
		clear sample_matrix;
		
		% Define the name of the subfolder
		data_subfolder_name = fullfile( Training_data_store_path, 'training_data_temp' );
		% Create the subfolder if it doesn't exist
		if ~exist( data_subfolder_name, 'dir' )
			mkdir( data_subfolder_name );
		end
		
		% Split the dataset into equal parts and save them to the subfolder
		for q = 1 : 1 : Num_replica
			training_data_index_q = ( q - 1 ) * 2 * N + 1 : q * 2 * N;
			training_data_part = data_set( training_data_index_q, : );
			% Save the current part to the subfolder
			save( fullfile( data_subfolder_name, sprintf( 'Training_data_part_%d.mat', q ) ), 'training_data_part' );
		end
		
		Langevin_sampling_test = 1;
		
	end
	