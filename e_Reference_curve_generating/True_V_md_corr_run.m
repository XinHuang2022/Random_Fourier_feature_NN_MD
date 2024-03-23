    
function Test_run = True_V_md_corr_run( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, Num_parallel_worker, L_x_vectorize, save_ref_curve, filename0, datadir0 )
	
	%{
    % First load the model trained with Tensorflow and Keras
	modelFolder_2 = './Model_save_HO_shallow_1';
	net = importTensorFlowNetwork( modelFolder_2, 'TargetNetwork', 'dlnetwork' );
	analyzeNetwork(net)
    %}
	
	% First load the model parameters
    %{
	omega = readmatrix( './parameter_set/K=256_paraset/omega_frequencies_K=256.txt' );
	eta = readmatrix( './parameter_set/K=256_paraset/weight_parameters_K=256.txt' );
	eta_re = eta( 1, : );
	eta_im = eta( 2, : );
    %}
	
		dirname = datadir0;
		% outfile = [ filename0,'.mat' ];

		corr_curve_info_1 = corr_function_curve( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, Num_parallel_worker, L_x_vectorize );
		correlation_md_record_mean_x = corr_curve_info_1( 1, : );
		correlation_md_record_std_x = corr_curve_info_1( 2, : );
		correlation_md_record_mean_x_2h = corr_curve_info_1( 3, : );
		
		correlation_md_record_mean_p = corr_curve_info_1( 4, : );
		correlation_md_record_std_p = corr_curve_info_1( 5, : );
		correlation_md_record_mean_p_2h = corr_curve_info_1( 6, : );
		
		% corr_curve_info_2 = corr_function_curve( 2 * h_x_sample, M_sp, tau, Delta_tau_plot );
		correlation_md_record_h_diff_x = abs( correlation_md_record_mean_x - correlation_md_record_mean_x_2h );
		correlation_md_record_h_diff_p = abs( correlation_md_record_mean_p - correlation_md_record_mean_p_2h );
		
		tau_values = ( 0 : Delta_tau_plot : tau );
		
		fig2 = figure( 2 );
		set( groot, 'defaultAxesTickLabelInterpreter', 'latex' ); 
		set( groot, 'defaultLegendInterpreter', 'latex' );
		plot( tau_values, correlation_md_record_mean_x );
		hold on 
		tau_values_2 = [ tau_values, fliplr( tau_values ) ];
		CI_upper_corr_curve = correlation_md_record_mean_x + correlation_md_record_std_x * 1.96 / sqrt( Num_replica );
		CI_lower_corr_curve = correlation_md_record_mean_x - correlation_md_record_std_x * 1.96 / sqrt( Num_replica );
		inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
		fill( tau_values_2, inBetween_statistic, 'cyan', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
		
		diff_upper_corr_curve = correlation_md_record_mean_x + correlation_md_record_h_diff_x;
		diff_lower_corr_curve = correlation_md_record_mean_x - correlation_md_record_h_diff_x;
		inBetween_h_diff = [ diff_lower_corr_curve, fliplr( diff_upper_corr_curve )];
		fill( tau_values_2, inBetween_h_diff, 'magenta', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
		
		% plot( tau_values, correlation_exact );
		% hold off
		% legend( 'MD with trained $\bar{V}(x)$', 'Exact value' );
		legend( 'MD with trained $\bar{V}(x)$', 'statistical CI', 'difference between h and 2h' );
		xlabel( 'correlation time $\tau$' )
		ylabel( '$\langle x_1(\tau),x_1(0)\rangle$-correlation' )
		
		
		fig3 = figure( 3 );
		set( groot, 'defaultAxesTickLabelInterpreter', 'latex' ); 
		set( groot, 'defaultLegendInterpreter', 'latex' );
		plot( tau_values, correlation_md_record_mean_p );
		hold on 
		tau_values_2 = [ tau_values, fliplr( tau_values ) ];
		CI_upper_corr_curve = correlation_md_record_mean_p + correlation_md_record_std_p * 1.96 / sqrt( Num_replica );
		CI_lower_corr_curve = correlation_md_record_mean_p - correlation_md_record_std_p * 1.96 / sqrt( Num_replica );
		inBetween_statistic = [ CI_lower_corr_curve, fliplr( CI_upper_corr_curve )];
		fill( tau_values_2, inBetween_statistic, 'cyan', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
		
		diff_upper_corr_curve = correlation_md_record_mean_p + correlation_md_record_h_diff_p;
		diff_lower_corr_curve = correlation_md_record_mean_p - correlation_md_record_h_diff_p;
		inBetween_h_diff = [ diff_lower_corr_curve, fliplr( diff_upper_corr_curve )];
		fill( tau_values_2, inBetween_h_diff, 'magenta', 'FaceAlpha', 0.2, 'LineStyle', 'none' );
		
		% plot( tau_values, correlation_exact );
		% hold off
		% legend( 'MD with trained $\bar{V}(x)$', 'Exact value' );
		legend( 'MD with trained $\bar{V}(x)$', 'statistical CI', 'difference between h and 2h' );
		xlabel( 'correlation time $\tau$' )
		ylabel( '$\langle p_1(\tau),p_1(0)\rangle$-correlation' )
		

		max_std_tau_curve = max( correlation_md_record_std_x );
		max_h_2h_diff_tau_curve = max( correlation_md_record_h_diff_x );

        max_std_tau_curve
        max_h_2h_diff_tau_curve
		%{
		diff_MD_exact = abs( correlation_exact - correlation_md_record );
		[ max_error, max_error_index ] = max( diff_MD_exact );
		fprintf( 'maximum error on x_1 correlation approximation is %4.3f, happening at tau = %3.2f\n', max_error, tau_values( max_error_index , 1 ) );
		
		abs_error = abs( correlation_md_record - correlation_exact );
		rel_error = abs( correlation_md_record - correlation_exact ) ./ correlation_exact;
		[max_rel_error, max_rel_error_id ] = max( rel_error );
		max_rel_error
		max_rel_error_id

		rel_error_tau_1 = abs_error( 1 + 1 / Delta_tau_plot, 1 ) / correlation_exact( 1 + 1 / Delta_tau_plot, 1 );
		rel_error_tau_1
		mean( abs_error )
		%}
		
		% Save the correlation curve data
		if( save_ref_curve )
			h_tau = h_x_sample;
			% save( 'correlation_md_curve_h_0.005.txt', 'correlation_md_record', '-ascii' );
			corr_curve_save_filename = sprintf( 'correlation_md_curve_h_%6.5f.mat', h_tau );
			% save('correlation_md_curve_h_0.001.mat', 'correlation_md_record');
			save( corr_curve_save_filename, 'corr_curve_info_1');
			
			folderName = 'Figures_save';
			if ~exist( folderName, 'dir' )
				mkdir( folderName );
			end
			corr_figure_save_filename_1 = sprintf( 'correlation_md_x_ref_h=%6.5f_M_sp=%d_N_rep=%d.fig', h_tau, M_sp, Num_replica );
			saveas( fig2, fullfile( folderName, corr_figure_save_filename_1 ) ); % Using figure number
			
			corr_figure_save_filename_2 = sprintf( 'correlation_md_p_ref_h=%6.5f_M_sp=%d_N_rep=%d.fig', h_tau, M_sp, Num_replica );
			saveas( fig3, fullfile( folderName, corr_figure_save_filename_2 ) ); % Using figure number
		
            fprintf( 1, '[END] PIMC_BB_sample_run completed: data saved to %s\n', corr_curve_save_filename );
        end
		
		
		
        Test_run = 1;
	end
	
	
	
	
	
	
	