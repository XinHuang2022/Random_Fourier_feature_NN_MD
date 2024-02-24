    
	clear all
	close all
	
	
	dim = 2;
	alpha = sqrt( 2 );
	gamma = 2;
	R_a = 2;
	R_c = 4;
	sigma_a = 1;
	sigma_c = 0.2;
	h_x_sample = 0.005;
	M_sp = 2^15;
	Num_replica = 32;
	Num_parallel_worker = 8;
	L_x_vectorize = 16;
	
	
	% Display prompt and read user input
	Mix_beta_true = input( 'Use single or mixed beta for sampling?\n For single beta please input 1, for mixed beta please input 2:\n ' );
	
	if( Mix_beta_true == 2 )
		beta_1 = input( 'Please enter the first beta value: \n' );
		fprintf( 'The first inverse temperature for sampling beta_1 = %3.2f\n', beta_1 );
		beta_2 = input( 'Please enter the second beta value: \n' );
		fprintf( 'The second inverse temperature for sampling beta_2 = %3.2f\n', beta_2 );
		fprintf( 'Sampling will be based on using two inverse temperature values.\n');
	elseif( Mix_beta_true == 1 )
		beta = input( 'Please enter the beta value for sampling: \n' );
		fprintf( 'The inverse temperature for sampling beta = %3.2f\n', beta );
		fprintf( 'Sampling will be based on using one inverse temperature value.\n');
	end
	
	if( Mix_beta_true == 2 )
		exitcode = Langevin_sampling_mix_beta( beta_1, beta_2, dim, alpha, gamma, R_a, R_c, sigma_a, sigma_c, h_x_sample, M_sp, Num_replica, Num_parallel_worker, L_x_vectorize )
	elseif( Mix_beta_true == 1 )
		exitcode = Langevin_sampling_one_beta( beta, dim, alpha, gamma, R_a, R_c, sigma_a, sigma_c, h_x_sample, M_sp, Num_replica, Num_parallel_worker, L_x_vectorize )
	end
	
	
	
	
	
	
	