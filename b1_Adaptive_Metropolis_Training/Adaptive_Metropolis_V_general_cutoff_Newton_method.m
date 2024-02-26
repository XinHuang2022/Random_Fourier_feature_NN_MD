   
    clear all
    close all
	
	dim = 2;        					% dim is the dimensionality
	N = 1 * 10^3;     					% N is the training data set size
	N_test = N;   					% N_test is the testing data set size
	lambda_1 = 0.01; 					% lambda is the Tikhonov regularization parameter
	lambda_2 = 0.001;
	
	M = 1;      					% M is the number of Metropolis iterations
	delta_step = ( 2.4^2 ) / dim;			
	delta = delta_step / 10;		% delta is the step length in Metropolis trial proposals
	gamma_expo = 3 * dim - 2; 				% gamma_expo is the exponent parameter in Metropolis test.
	delta_lr_omega = 0.0001;				% delta_lr is the learning rate for SGD with training on the frequencies
    % delta_lr_omega = 0.01;
	
	% K_values = [ 8, 16, 32, 64, 128, 256, 512, 1024 ];
	% K_values = [ 8, 16, 32, 64, 128 ];
    % K_values = [ 8, 16, 32, 64, 128, 256, 512 ];
    K_values = 64;

	num_K_values = length( K_values );
	Q = 2;       					% Q is the number of replica for evaluating the error bar of the loss function
	
	Errors_K_rec_1 = zeros( num_K_values, 2 );    % Errors_K_rec_1 stores the mean and error bars of loss in the adaptive Metropolis algorithm
	Errors_K_rec_2 = zeros( num_K_values, 2 );    % Errors_K_rec_2 stores the mean and error bars of loss in the fixed random feature algorithm
	
	Error_rep_q_1 = zeros( num_K_values, Q );
	Error_rep_q_2 = zeros( num_K_values, Q );
	
	R_a = 2;					% Specify the parameters for the smoothed cut-off function \chi_1( r )
	% R_b = 3;
	R_c = 4;
	sigma_a = 1;
	sigma_c = 0.2;
	
	Mix_beta_flag = input( 'Use single or mixed beta for training data?\n For single beta please input 1, for mixed beta please input 2:\n ' );
	Use_laptop_flag = input( 'Using personal laptop or the office desktop?\n For personal laptop please input 1, for office desktop please input 2:\n ' );
	M_sp = 2^15;    % M_sp is the sample set size in the data generation pipeline with Langevin dynamics sampling
	Num_replica = 32;    % Num_replica is the number of independent training data set from the data generation pipeline
	Delta_t_sample = 2;    % Delta_t_sample is the time difference between two samples in the Langevin dynamics sampling step
	h_x_sample = 0.005;    % h_x_sample is the time step size in the Langevin dynamics sampling 
	
	
	%{
	if( Mix_beta_flag == 2 )
		if( Use_laptop_flag == 1 )
			csv_file_path_1 = 'G:\Research\Projects\Path_Integral_NN_MD\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_mix_betadim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
		elseif( Use_laptop_flag == 2 )
			csv_file_path_1 = 'C:\Users\ander\Documents\MATLAB\Path_Integral_MD_NN\24021901\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_mix_betadim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
		end
	elseif( Mix_beta_flag == 1 )
		if( Use_laptop_flag == 1 )
			csv_file_path_1 = 'G:\Research\Projects\Path_Integral_NN_MD\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_one_beta_dim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
		elseif( Use_laptop_flag == 2 )
			csv_file_path_1 = 'C:\Users\ander\Documents\MATLAB\Path_Integral_MD_NN\24021901\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_one_beta_dim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
		end
	end
	%}
	if( Mix_beta_flag == 2 )
		csv_file_path_1 = "..\d1_Visualisation_and_Data_Storage\Data_saved\Sample_store\Langevin_sample_store_V_general_mix_beta" + "_dim=" + dim + "_Msp=" + M_sp + "_Num_rep=" + Num_replica + "_Delta_t_sample=" + Delta_t_sample + "_dt=" + h_x_sample + ".csv"
	elseif( Mix_beta_flag == 1 )
		csv_file_path_1 = "..\d1_Visualisation_and_Data_Storage\Data_saved\Sample_store\Langevin_sample_store_V_general_one_beta" + "_dim=" + dim + "_Msp=" + M_sp + "_Num_rep=" + Num_replica + "_Delta_t_sample=" + Delta_t_sample + "_dt=" + h_x_sample + ".csv"
	end
	% csv_file_path_1 = 'G:\Research\Projects\Path_Integral_NN_MD\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_mix_betadim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
	% csv_file_path_2 = 'G:\Research\Projects\Path_Integral_NN_MD\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_dim=2_Msp=1000000_Delta_t_sample=1_dt=0.01_beta=1.csv';
	% csv_file_path_1 = 'C:\Users\ander\Documents\Python_space\Projects\potential_training_6\tests\24012902\2_2_V_general\b2_SGD_Training\Langevin_sample_store_V_general_dim=2_Msp=1000000_Delta_t_sample=1_dt=0.01_beta=0.2.csv';
	% csv_file_path_2 = 'C:\Users\ander\Documents\Python_space\Projects\potential_training_6\tests\24012902\2_2_V_general\b2_SGD_Training\Langevin_sample_store_V_general_dim=2_Msp=1000000_Delta_t_sample=1_dt=0.01_beta=1.csv';
	% csv_file_path_1 = '/cfs/klemming/home/x/xinhuang/Private/2024021801/Langevin_sample_store_V_general_mix_betadim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
	% csv_file_path_2 = '/cfs/klemming/home/x/xinhuang/Private/2024020101/Langevin_sample_store_V_general_dim=2_Msp=1000000_Delta_t_sample=1_dt=0.01_beta=1.csv';
	% csv_file_path_1 = 'C:\Users\ander\Documents\MATLAB\Path_Integral_MD_NN\24021901\24020801\2_2_V_general\a_Data_sampling\Sample_store\Langevin_sample_store_V_general_mix_betadim=2_Msp=32768_Num_rep=32_Delta_t_sample=2_dt=0.005.csv';
    data_set_1 = readmatrix( csv_file_path_1 );
	% data_set_2 = readmatrix( csv_file_path_2 );
	
	% Reshuffle on the data set with mixed beta values
	% data_set_concat = [ data_set_1( 1 : end, : ); data_set_2( 1 : end, : ) ];
	% [ N_row, N_col ] = size( data_set_concat );
	[ N_row, N_col ] = size( data_set_1 );
	rand_indices = randperm( N_row );
	data_set = zeros( N_row, N_col );
	data_set( rand_indices, : )  = data_set_1;
	
    % parpool( 8 );
    tic 
	for q = 1 : 1 : Q
    % parfor q = 1 : Q
		
		training_data_index_q = ( q - 1 ) * 2 * N + 1 : ( q - 1 ) * 2 * N + N;
		x_data = data_set( training_data_index_q, 1 : 2 );
		v_data = data_set( training_data_index_q, 3 );
		v_prime_1_data = data_set( training_data_index_q, 4 );
		v_prime_2_data = data_set( training_data_index_q, 5 );
		y_data = [ v_data; v_prime_1_data; v_prime_2_data ];
		
		testing_data_index_q = ( q - 1 ) * 2 * N + N + 1 : q * 2 * N;
		x_data_test = data_set( testing_data_index_q, 1 : 2 );
		v_data_test = data_set( testing_data_index_q, 3 );
		v_prime_1_data_test = data_set( testing_data_index_q, 4 );
		v_prime_2_data_test = data_set( testing_data_index_q, 5 );
		y_data_test = [ v_data_test; v_prime_1_data_test; v_prime_2_data_test ];
		
		N_num_Newton_beta = 50;
		
		for j = 1 : 1 : num_K_values
			
			K = K_values( 1, j );
            rng( 126 )
            omega_AM_0 = 3 * randn( K, dim );
            omega_GD_0 = 3 * randn( K, dim );
			
			omega = omega_AM_0;
			S_v_mat = exp( 1i * ( x_data * omega' ) );
			omega_dim_1 = 1i * omega( :, 1 );
			S_v_grad_1_mat = ( omega_dim_1.' ) .* S_v_mat;
			omega_dim_2 = 1i * omega( :, 2 );
			S_v_grad_2_mat = ( omega_dim_2.' ) .* S_v_mat;
			S_mat = [ S_v_mat; S_v_grad_1_mat; S_v_grad_2_mat ];
			beta_hat = ( S_mat' * S_mat + lambda_1 * N * diag( ones( K, 1 ) ) ) \ ( S_mat' * y_data );
			for m = 1 : 1 : M
				
				r_normal_step = randn( K, dim );
				omega_prime = omega + delta * r_normal_step;
				S_v_mat_prime = exp( 1i * ( x_data * omega_prime' ) );
				omega_prime_dim_1 = 1i * omega_prime( :, 1 );
				S_v_grad_1_mat_prime = ( omega_prime_dim_1.' ) .* S_v_mat_prime;
				omega_prime_dim_2 = 1i * omega_prime( :, 2 );
				S_v_grad_2_mat_prime = ( omega_prime_dim_2.' ) .* S_v_mat_prime;
				S_mat_prime = [ S_v_mat_prime; S_v_grad_1_mat_prime; S_v_grad_2_mat_prime ];
				
				% beta_hat_prime = ( S_mat_prime' * S_mat_prime + lambda * N * diag( ones( K, 1 ) ) ) \ ( S_mat_prime' * y_data );
				beta_hat_prime = Newton_Raphson_beta( S_mat_prime, beta_hat, y_data, N, K, lambda_1, lambda_2, N_num_Newton_beta );
				
				r_uniform = rand( K, 1 );
				beta_k_increase = r_uniform < ( ( beta_hat_prime ./ beta_hat ).^gamma_expo );
				omega_new = beta_k_increase .* omega_prime + ( 1 - beta_k_increase ) .* omega;
				omega = omega_new;
				beta_hat = beta_hat_prime;
			
			end
			
			S_v_mat_final = exp( 1i * ( x_data * omega' ) );
			omega_final_dim_1 = 1i * omega( :, 1 );
			S_v_grad_1_mat_final = ( omega_final_dim_1.' ) .* S_v_mat_final;
			omega_final_dim_2 = 1i * omega( :, 2 );
			S_v_grad_2_mat_final = ( omega_final_dim_2.' ) .* S_v_mat_final;
			S_mat_final = [ S_v_mat_final; S_v_grad_1_mat_final; S_v_grad_2_mat_final ];
			
			% beta_hat_final = ( S_mat_final' * S_mat_final + lambda * N * diag( ones( K, 1 ) ) ) \ ( S_mat_final' * y_data );
			beta_hat_final = Newton_Raphson_beta( S_mat_final, beta_hat, y_data, N, K, lambda_1, lambda_2, N_num_Newton_beta );
			
			S_v_mat_test = exp( 1i * ( x_data_test * omega' ) );
			S_v_grad_1_mat_test = ( omega_final_dim_1.' ) .* S_v_mat_test;
			S_v_grad_2_mat_test = ( omega_final_dim_2.' ) .* S_v_mat_test;
			S_mat_test = [ S_v_mat_test; S_v_grad_1_mat_test; S_v_grad_2_mat_test ];
			
			% loss_K_adaptive = sum( abs( S_mat_test * beta_hat_final - y_data_test ).^2, 1 ) + lambda_1 * norm( beta_hat_final, 2 )^2;
			% loss_K_adaptive = sum( abs( S_mat_test * beta_hat_final - y_data_test ).^2, 1 ) + lambda_1 * norm( beta_hat_final, 2 )^2 + lambda_2 * ( norm( beta_hat_final, 2 )^2 )^2;
			loss_K_adaptive = sum( abs( S_mat_test * beta_hat_final - y_data_test ).^2, 1 );
			Error_rep_q_1( j, q ) = loss_K_adaptive / N_test;
			
			if ~exist( 'frequencies', 'dir' )			% Create 'frequencies' folder if it doesn't exist
				mkdir('frequencies');
			end
			if ~exist( 'amplitudes', 'dir' )			% Create 'amplitudes' folder if it doesn't exist
				mkdir('amplitudes');
			end
			
			filename_omega_AM = sprintf( 'frequency_omega_AM_K=%d_q=%d.csv', K, q );
			filepath_omega_AM = fullfile( 'frequencies', filename_omega_AM );  % Full path including 'frequencies' subfolder
			writematrix( omega', filepath_omega_AM );
			beta_AM_real = real( beta_hat_final );
			beta_AM_imag = imag( beta_hat_final );
			beta_AM_final = [ beta_AM_real, beta_AM_imag ];
			filename_eta_AM = sprintf( 'amplitude_eta_AM_K=%d_q=%d.csv', K, q );
			filepath_eta_AM = fullfile( 'amplitudes', filename_eta_AM );  % Full path including 'amplitudes' subfolder
			writematrix( beta_AM_final', filepath_eta_AM );
			
			
			
			omega_gd = omega_GD_0;
			% S_mat_sgd = exp( 1i * ( x_data * omega_sgd' ) );
			S_v_mat_gd = exp( 1i * ( x_data * omega_gd' ) );
			omega_dim_1_gd = 1i * omega_gd( :, 1 );
			S_v_grad_1_mat_gd = ( omega_dim_1_gd.' ) .* S_v_mat_gd;
			omega_dim_2_gd = 1i * omega_gd( :, 2 );
			S_v_grad_2_mat_gd = ( omega_dim_2_gd.' ) .* S_v_mat_gd;
			S_mat_gd = [ S_v_mat_gd; S_v_grad_1_mat_gd; S_v_grad_2_mat_gd ];
			% beta_hat_sgd = ( S_mat_sgd' * S_mat_sgd + lambda * N * diag( ones( K, 1 ) ) ) \ ( S_mat_sgd' * y_data );
			beta_hat_gd_init = ones( K, 1 );
			beta_hat_gd_0 = Newton_Raphson_beta( S_mat_gd, beta_hat_gd_init, y_data, N, K, lambda_1, lambda_2, N_num_Newton_beta );
			beta_hat_gd = beta_hat_gd_0;
            beta_hat_gd
            for m = 1 : 1 : M
			    
                %{
				dL_dbeta = ( 1 / N ) * 2 * S_mat_sgd' * ( S_mat_sgd * beta_hat_sgd - y_data ) + 2 * lambda * beta_hat_sgd;
				dL_domega = ( 1 / N ) * 2 * ( S_mat_sgd' * ( ( S_mat_sgd * beta_hat_sgd - y_data ) .* ( 1i * x_data ) ) ) .* beta_hat_sgd;
				
				beta_hat_sgd_new = beta_hat_sgd - delta_lr * dL_dbeta;
				omega_sgd_new = omega_sgd - delta_lr * dL_domega;
				
				beta_hat_sgd = beta_hat_sgd_new;
				omega_sgd = omega_sgd_new;
                %}
                % dL_domega_1 = ( 1 / N ) * 2 * real( ( ( S_mat_sgd * beta_hat_sgd - y_data ).' * ( conj( S_mat_sgd ) .* ( -1i * x_data( :, 1 ) ) ) ) .* ( beta_hat_sgd.' ) );
				% dL_domega_2 = ( 1 / N ) * 2 * real( ( ( S_mat_sgd * beta_hat_sgd - y_data ).' * ( conj( S_mat_sgd ) .* ( -1i * x_data( :, 2 ) ) ) ) .* ( beta_hat_sgd.' ) );
				
				dL1_domega_1 = ( 2 / N ) * real( ( ( S_v_mat_gd * beta_hat_gd - v_data ).' * ( conj( S_v_mat_gd ) .* ( -1i * x_data( :, 1 ) ) ) ) .* ( beta_hat_gd.' ) );
				dL1_domega_2 = ( 2 / N ) * real( ( ( S_v_mat_gd * beta_hat_gd - v_data ).' * ( conj( S_v_mat_gd ) .* ( -1i * x_data( :, 2 ) ) ) ) .* ( beta_hat_gd.' ) );
				dL21_domega_1 = ( 2 / N ) * real( ( ( -1i ) * ( conj( beta_hat_gd ) .* S_v_mat_gd' ) + ( omega_gd( :, 1 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 1 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_1_mat_gd * beta_hat_gd - v_prime_1_data ) );
                dL21_domega_2 = ( 2 / N ) * real( ( ( omega_gd( :, 1 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 2 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_1_mat_gd * beta_hat_gd - v_prime_1_data ) );
                dL22_domega_1 = ( 2 / N ) * real( ( ( omega_gd( :, 2 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 1 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_2_mat_gd * beta_hat_gd - v_prime_2_data ) );
                dL22_domega_2 = ( 2 / N ) * real( ( ( -1i ) * ( conj( beta_hat_gd ) .* S_v_mat_gd' ) + ( omega_gd( :, 2 ) .* conj( beta_hat_gd ) ) .* ( -x_data( :, 2 ).' .* S_v_mat_gd' ) ) * ( S_v_grad_2_mat_gd * beta_hat_gd - v_prime_2_data ) );
                
				dL_domega_1 = dL1_domega_1.' + dL21_domega_1 + dL22_domega_1;
				dL_domega_2 = dL1_domega_2.' + dL21_domega_2 + dL22_domega_2;
				
                dL_domega = [ dL_domega_1, dL_domega_2 ];
                max( max( abs( dL_domega ) ) )

				% beta_hat_sgd_new = beta_hat_sgd - delta_lr_beta * dL_dbeta;
				omega_gd_new = omega_gd - delta_lr_omega * dL_domega;
				
				% beta_hat_sgd = beta_hat_sgd_new;
				omega_gd = omega_gd_new;
				% S_mat_sgd = exp( 1i * ( x_data * omega_sgd' ) );
				S_v_mat_gd = exp( 1i * ( x_data * omega_gd' ) );
				omega_dim_1_gd = 1i * omega_gd( :, 1 );
				S_v_grad_1_mat_gd = ( omega_dim_1_gd.' ) .* S_v_mat_gd;
				omega_dim_2_gd = 1i * omega_gd( :, 2 );
				S_v_grad_2_mat_gd = ( omega_dim_2_gd.' ) .* S_v_mat_gd;
				S_mat_gd = [ S_v_mat_gd; S_v_grad_1_mat_gd; S_v_grad_2_mat_gd ];
			    % beta_hat_sgd = ( S_mat_sgd' * S_mat_sgd + lambda * N * diag( ones( K, 1 ) ) ) \ ( S_mat_sgd' * y_data );
				beta_hat_gd = Newton_Raphson_beta( S_mat_gd, beta_hat_gd_0, y_data, N, K, lambda_1, lambda_2, N_num_Newton_beta );
				beta_hat_gd_0 = beta_hat_gd;
				
			end
				
			% S_mat_gd = exp( 1i * ( x_data * omega_sgd' ) );
			% beta_hat_sgd = ( S_mat_sgd' * S_mat_sgd + lambda * N * diag( ones( K, 1 ) ) ) \ ( S_mat_sgd' * y_data );
			% beta_hat_sgd = Newton_Raphson_beta( S_mat_sgd, beta_hat_sgd_0, y_data, N, K, lambda_1, lambda_2, N_num_Newton_beta );
			% S_mat_gd_test = exp( 1i * ( x_data_test * omega_gd' ) );

            S_v_mat_gd_test = exp( 1i * ( x_data_test * omega_gd' ) );
			S_v_grad_1_mat_gd_test = ( omega_dim_1_gd.' ) .* S_v_mat_gd_test;
			S_v_grad_2_mat_gd_test = ( omega_dim_2_gd.' ) .* S_v_mat_gd_test;
			S_mat_gd_test = [ S_v_mat_gd_test; S_v_grad_1_mat_gd_test; S_v_grad_2_mat_gd_test ];

			% loss_K_gd = sum( abs( S_mat_gd_test * beta_hat_gd - y_data_test ).^2, 1 ) + lambda_1 * norm( beta_hat_final, 2 )^2;
			% loss_K_gd = sum( abs( S_mat_gd_test * beta_hat_gd - y_data_test ).^2, 1 ) + lambda_1 * norm( beta_hat_gd, 2 )^2 + lambda_2 * ( norm( beta_hat_gd, 2 )^2 )^2;
			loss_K_gd = sum( abs( S_mat_gd_test * beta_hat_gd - y_data_test ).^2, 1 );
			Error_rep_q_2( j, q ) = loss_K_gd / N_test;
			
			filename_omega_GD = sprintf( 'frequency_omega_GD_K=%d_q=%d.csv', K, q );
			filepath_omega_GD = fullfile( 'frequencies', filename_omega_GD );  % Full path including 'frequencies' subfolder
			writematrix( omega_gd', filepath_omega_GD );
			beta_GD_real = real( beta_hat_gd );
			beta_GD_imag = imag( beta_hat_gd );
			beta_GD_final = [ beta_GD_real, beta_GD_imag ];
			filename_eta_GD = sprintf( 'amplitude_eta_GD_K=%d_q=%d.csv', K, q );
			filepath_eta_GD = fullfile( 'amplitudes', filename_eta_GD );  % Full path including 'amplitudes' subfolder
			writematrix( beta_GD_final', filepath_eta_GD );
			
			% fprintf( 'K = %d completed.\n', K );
        end

        fprintf( '%d iterations for q.\n', q );
		
    end
    poolobj = gcp( 'nocreate' );
    delete( poolobj );

    toc
	
	Errors_K_rec_1( :, 1 ) = mean( Error_rep_q_1, 2 );
	Errors_K_rec_1( :, 2 ) = std( Error_rep_q_1, 0, 2 );
	
	Errors_K_rec_2( :, 1 ) = mean( Error_rep_q_2, 2 );
	Errors_K_rec_2( :, 2 ) = std( Error_rep_q_2, 0, 2 );

    K_values_sgd = [ 64, 128, 256, 512, 1024 ];
    Errors_K_rec_3 = [ 0.4995, 0.1517, 0.0910, 0.0556, 0.0309 ];
	
	figure( 1 )
	errorbar( K_values, Errors_K_rec_1( :, 1 ), Errors_K_rec_1( :, 2 ) * 1.96 / sqrt( Q ), '-*' );
    % loglog( K_values, Errors_K_rec_1( :, 1 ), '-*' );
	hold on
	errorbar( K_values, Errors_K_rec_2( :, 1 ), Errors_K_rec_2( :, 2 ) * 1.96 / sqrt( Q ), '-o' );
    % loglog( K_values, Errors_K_rec_2( :, 1 ), '-o' );
    % loglog( K_values_sgd, Errors_K_rec_3, '-+' )
	loglog( K_values, 25 * K_values.^(-1) );
	hold off
	% legend( 'loss for adaptive Metropolis algorithm', 'loss for Gradient Descent on frequency', 'loss for SGD on frequency and amplitude', 'reference $\mathcal{O}(K^{-1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
	legend( 'loss for adaptive Metropolis algorithm', 'loss for Gradient Descent on frequency', 'reference $\mathcal{O}(K^{-1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
    % legend( 'loss for adaptive Metropolis algorithm', 'reference $\mathcal{O}(K^{-1})$', 'fontsize', 20, 'interpreter', 'latex' ); 
	xlabel( '$K$', 'fontsize', 20, 'interpreter', 'latex' );
    ylabel( 'Loss on the test set', 'fontsize', 20, 'interpreter', 'latex' );
    % title( 'Log-log plot of testing loss, target function $f(x)=e^{-\frac{|x|^{2}}{2}}\,\mathrm{Si}(\frac{x}{a})$, $a=1$, $J=10^3$, $d=1$', 'fontsize', 20, 'interpreter', 'latex' )
	% title( 'Log-log plot of testing loss, target function $v_1(x)=\frac{|x_1|^2}{2}\,\chi_1(|x|)$, $J=10^3$, $d=2$', 'fontsize', 20, 'interpreter', 'latex' )
    title_string = sprintf( 'Log-log plot of testing loss, more complicated target function, $J=%d$, $d=%d$, $Q=%d$', N, dim, Q );
    title( title_string, 'fontsize', 20, 'interpreter', 'latex' )
    set( gca, 'YScale', 'log');
    set( gca, 'XScale', 'log');
	
    
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
    
		
		
		
		
	
	
	
	
	
	