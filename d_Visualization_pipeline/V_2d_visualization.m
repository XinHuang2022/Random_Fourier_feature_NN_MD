    
	function V_2d_plot_status = V_2d_visualization( R_a, sigma_a, R_b, sigma_b, K, Num_replica, Training_pipeline_path )
		
		x = -5 : 0.02 : 5;
		y = -5 : 0.02 : 5;
		[ X, Y ] = meshgrid( x, y ); 
		
		% Compute the function values for the grid of (x, y) points
		Z_true = V_true_2d_plot( X, Y, R_a, sigma_a );
		Z_approx = V_approx_2d_plot( X, Y, R_a, R_b, sigma_a, sigma_b, K, Num_replica, Training_pipeline_path );

		% Create a 3D surface plot
		fig9 = figure( 9 );
		surf( X, Y, Z_true, 'EdgeColor', 'none' );
		% contour( X, Y, Z_true, 'EdgeColor', 'none' );

		xlabel( 'x_1-axis' );
		ylabel( 'x_2-axis' );
		zlabel( 'V(x_1, x_2)' );
		title( '3D Plot of the Function V = 0.5 x_1^2 + 0.5 \alpha x_2^2 + \gamma sin( x_1 x_2 ) * \chi_{-1}(r)');

		colormap( 'jet' );
		colorbar;       
		view( 0, 90 ); 
		

		fig10 = figure( 10 );
		surf( X, Y, Z_approx, 'EdgeColor', 'none' );

		xlabel( 'x_1-axis' );
		ylabel( 'x_2-axis' );
		zlabel( '\bar{V}(x_1, x_2)' );
		title( '3D Plot of the approximated function \bar{V}');

		colormap( 'jet' );
		colorbar;        
		view( 0, 90 ); 


		fig11 = figure( 11 );
		surf( X, Y, abs( Z_approx - Z_true ), 'EdgeColor', 'none' );

		xlabel( 'x_1-axis' );
		ylabel( 'x_2-axis' );
		zlabel( '|V(x_1, x_2) - \bar{V}(x_1, x_2)|' );
		title( '3D Plot of the difference between V(x) and \bar{V}(x)');

		colormap( 'jet' ); 
		colorbar;       
		view( 0, 90 ); 
		
		folderName = fullfile( 'd_Visualization_pipeline', 'V_figures_save' );
		if ~exist( folderName, 'dir' )
			mkdir( folderName );
		end
		
		V_true_plot_save_filename = sprintf( 'V_true_plot_2d.fig' );
		saveas( fig9, fullfile( folderName, V_true_plot_save_filename ) );
		
		V_reconstruct_plot_save_filename = sprintf( 'V_reconstruct_AM_K=%d_Nrep=%d.fig', K, Num_replica );
		saveas( fig10, fullfile( folderName, V_reconstruct_plot_save_filename ) ); 
		
		V_difference_plot_save_filename = sprintf( 'V_pointwise_diff_K=%d_Nrep=%d.fig', K, Num_replica );
		saveas( fig11, fullfile( folderName, V_difference_plot_save_filename ) ); 
		
		V_2d_plot_status = 1;
		
	end
	