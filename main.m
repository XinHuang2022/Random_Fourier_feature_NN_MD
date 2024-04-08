	
	clear all
	close all
	
	% Add the path for each pipeline
	addpath( 'a_Data_sampling_pipeline' );
	addpath( 'b_Training_pipeline' );
	addpath( 'c_Inference_pipeline' );
	addpath( 'd_Visualization_pipeline' );
	addpath( 'e_Reference_curve_generating' );
	
	
	%%  Data sampling and collection pipeline
	dim = 2;
	alpha = sqrt( 2 );    gamma = 2;   					% alpha and gamma are the parameters in the potential function 
	R_a = 2;  R_b = 3;  R_c = 4;						% R_a, R_b, and R_c are the radius parameters corresponding to the cut-off functions \chi_{-1}( r ), \chi_0( r ), and \chi_1( r )
	sigma_a = 1;  sigma_b = 0.2;  sigma_c = 0.2;				% sigma_a, sigma_b, and sigma_c are the sharpness parameters corresponding to the cut-off functions \chi_{-1}( r ), \chi_0( r ), and \chi_1( r )
	h_x_sample = 0.005;							% h_x_sample is the time step length for the Langevin dynamics sampling
	J = 10^5;								% J is the training data set size for each independent realization of training and inference
    	Num_replica = 32;							% Num_replica is the number of independent replicas to evaluate the mean with respect to independent data set and frequency parameter sets
	
	Total_sample_power_index = ceil( log( J * Num_replica ) / log( 2 ) ) + 1;
	
	M_sp = 2^Total_sample_power_index;  					% M_sp is the total sample size for Gibbs density with Langevin dynamics sampling
	Num_parallel_worker = 128;						% Num_parallel_worker is the number of available CPU cores for the parallel pool            
	
	Mix_beta_true = 2;							% Mix_beta_true = 2 to use two inverse temperatures for data sampling, and Mix_beta_true = 1 to use only one inverse temperature
	if( Mix_beta_true == 2 )
		beta_1 = 1;
		beta_2 = 0.3;
	elseif( Mix_beta_true == 1 )
		beta = 1;
	end
	
    	current_dir = pwd;
    	Training_pipeline_path = fullfile( current_dir, 'b_Training_pipeline' );
	if( Mix_beta_true == 2 )
		Data_sampling_Flag = Langevin_sampling_mix_beta( beta_1, beta_2, dim, alpha, gamma, R_a, R_c, sigma_a, sigma_c, h_x_sample, M_sp, Num_replica, Num_parallel_worker, J, Training_pipeline_path );
	elseif( Mix_beta_true == 1 )
		Data_sampling_Flag = Langevin_sampling_one_beta( beta, dim, alpha, gamma, R_a, R_c, sigma_a, sigma_c, h_x_sample, M_sp, Num_replica, Num_parallel_worker, J, Training_pipeline_path );
	end
	
	%%  Training pipeline
	lambda_1 = 0.01; 							% lambda_1 is the Tikhonov regularization parameter
	lambda_2 = 0.001;							% lambda_2 is the penalty parameter corresponding to the square of sum of squared amplitudes
	lambda_3 = 0.01;							% lambda_3 is the parameter for the last penalty with restrictions on the norm of amplitudes exceeding a constant C_const_bound
	C_const_bound = 100;
	
    	K_values = [ 16, 32, 64, 128, 256, 512, 1024 ];
	
	Training_Flag = Training_Adaptive_Metropolis( dim, J, lambda_1, lambda_2, lambda_3, C_const_bound, K_values, Num_replica, Training_pipeline_path );
	
	%%  Correlation function inference pipeline
	M_MC = 2^16;								% M_MC is the number of samples for the Monte Carlo integral to evaluate the approximated correlation function using the trained potential function
	tau = 2;								% tau is the largest correlation time when studying the correlation function and plotting the correlation curve
	Delta_tau_plot = 0.1;							% Delta_tau_plot is the distance between two points in the plot of the correlation function curve

	Inference_Flag = Test_approx_V_md_corr( h_x_sample, M_MC, Num_replica, tau, Delta_tau_plot, Num_parallel_worker, R_a, sigma_a, R_b, sigma_b, K_values, dim, alpha, gamma, J, Training_pipeline_path );
	
	%%  Visualization pipeline
	Visualization_Flag = V_2d_visualization( R_a, sigma_a, R_b, sigma_b, max( K_values ), Num_replica, Training_pipeline_path );
	
	
	
	
	
	
