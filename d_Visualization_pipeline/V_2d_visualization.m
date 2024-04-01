    
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
        %{
		contourf( X, Y, Z_true, 20, 'LineStyle', 'none' ); % Filled contour plot
		hold on;
		contour( X, Y, Z_true, 5, 'LineColor', 'k' ); % Contour lines
		colorbar;
		%}
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
		zlabel( '\bar{v}_r(x_1, x_2)' );
		title( '3D Plot of the approximated function \bar{V}');

		colormap( 'jet' );
		colorbar;        
		view( 0, 90 ); 
		% view( -37.5, 30 )


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
		
		
        % Generating the new figure with histogram of the frequency parameters
	    omega_store_AM = zeros( 2, K * Num_replica );
	    omega_store_GD = zeros( 2, K * Num_replica );
	    for q = 1 : 1 : Num_replica
            filename_omega_AM = sprintf( 'frequency_omega_AM_K=%d_q=%d.csv', K, q );
		    filepath_omega_AM = fullfile( Training_pipeline_path, 'frequencies', filename_omega_AM );  % Full path including 'frequencies' subfolder
		    omega_AM = readmatrix( filepath_omega_AM );
		    omega_store_AM( :, ( q - 1 ) * K + 1 : q * K ) = omega_AM;
		    
		    filename_omega_GD = sprintf( 'frequency_omega_GD_K=%d_q=%d.csv', K, q );
		    filepath_omega_GD = fullfile( Training_pipeline_path, 'frequencies', filename_omega_GD );  % Full path including 'frequencies' subfolder
		    omega_GD = readmatrix( filepath_omega_GD );
		    omega_store_GD( :, ( q - 1 ) * K + 1 : q * K ) = omega_GD;
	    end
	    
        % Plot of the empirical density of frequencies by AM
	    fig16 = figure( 16 );
	    nbins = [ 200, 200 ];
	    % Convert histogram counts to density values using histcounts2
	    [ counts, xedge, yedge ] = histcounts2( omega_store_AM( 1, : ), omega_store_AM( 2, : ), nbins, 'Normalization', 'pdf');
	    % Plot 3D histogram using surf
        xedge_new = xedge( 1 : end - 1 ) + ( xedge( 1, 2 ) - xedge( 1, 1 ) ) / 2;
        yedge_new = yedge( 1 : end - 1 ) + ( yedge( 1, 2 ) - yedge( 1, 1 ) ) / 2;
        [ Omega_1, Omega_2 ] = meshgrid( xedge_new, yedge_new ); 
        % counts = counts / ( ( xedge( 1, 2 ) - xedge( 1, 1 ) ) * ( yedge( 1, 2 ) - yedge( 1, 1 ) ) );
	    surf( Omega_1, Omega_2, counts, 'EdgeColor', 'none' );
    
	    xlabel('\omega_1')
	    ylabel('\omega_2')
	    zlabel('Density');
	    title('3D Histogram of empirical density function by AM ');
		fig_width = 5.5; % Width of the figure in inches
		fig_height = 4; % Height of the figure in inches
        set( gcf, 'Units', 'inches', 'Position', [ 0, 0, fig_width, fig_height ] ) ;
	    colorbar; 
	    view( 0, 90 );
        
         % Plot of the empirical density of frequencies by GD
        fig17 = figure( 17 );
	    nbins = [ 200, 200 ];
	    [ counts, xedge, yedge ] = histcounts2( omega_store_GD( 1, : ), omega_store_GD( 2, : ), nbins, 'Normalization', 'pdf' );
        xedge_new = xedge( 1 : end - 1 ) + ( xedge( 1, 2 ) - xedge( 1, 1 ) ) / 2;
        yedge_new = yedge( 1 : end - 1 ) + ( yedge( 1, 2 ) - yedge( 1, 1 ) ) / 2;
        [ Omega_1, Omega_2 ] = meshgrid( xedge_new, yedge_new ); 
	    surf( Omega_1, Omega_2, counts, 'EdgeColor', 'none' );
    
	    xlabel('\omega_1')
	    ylabel('\omega_2')
	    zlabel('Density');
	    title('3D Histogram of empirical density function by GD');
		set( gcf, 'Units', 'inches', 'Position', [ 0, 0, fig_width, fig_height ] ) ;
	    colorbar; 
	    view( 0, 90 );
	    
	    optimal_p_vals = zeros( 200, 200 );
	    h_integral = 0.1;
	    x1_points = -5 : h_integral : 5;
	    x2_points = -5 : h_integral : 5;
	    N_quad = length( x1_points );
	    weight_mat = ones( N_quad, N_quad );
	    for i_c = 2 : 1 : N_quad - 1
		    weight_mat( 1, i_c ) = 0.5;
		    weight_mat( N_quad, i_c ) = 0.5;
		    weight_mat( i_c, 1 ) = 0.5;
		    weight_mat( i_c, N_quad ) = 0.5;
	    end
	    weight_mat( 1, 1 ) = 0.25;
	    weight_mat( 1, N_quad ) = 0.25;
	    weight_mat( N_quad, 1 ) = 0.25;
	    weight_mat( N_quad, N_quad ) = 0.25;
	    weight_mat = weight_mat * h_integral^2;
	    
	    [ X1_grid_integ, X2_grid_integ ] = meshgrid( x1_points, x2_points );
	    for i_a = 1 : 1 : 200
		    omega_1 = xedge_new( 1, i_a );
		    for i_b = 1 : 1 : 200
			    omega_2 = yedge_new( 1, i_b );
			    v1_FT_mesh = v1_FT_func( X1_grid_integ, X2_grid_integ, omega_1, omega_2 );
			    optimal_p_vals( i_a, i_b ) = abs( sum( weight_mat .* v1_FT_mesh, "all" ) ) * sqrt( 1 + omega_1^2 + omega_2^2 );
		    end
	    end
	    
	    optimal_p_vals_sum = sum( optimal_p_vals, "all" ) * ( xedge( 1, 2 ) - xedge( 1, 1 ) ) * ( yedge( 1, 2 ) - yedge( 1, 1 ) );
	    optimal_p_vals = optimal_p_vals / optimal_p_vals_sum;
	    
	    fig18 = figure( 18 );
	    surf( Omega_1, Omega_2, optimal_p_vals, 'EdgeColor', 'none' );
	    xlabel('\omega_1')
	    ylabel('\omega_2')
	    zlabel('Density');
	    title('3D Histogram of optimal density function');
		set( gcf, 'Units', 'inches', 'Position', [ 0, 0, fig_width, fig_height ] ) ;
        colorbar; 
        view( 0, 90 )
		
		p_omega_emp_AM_filename = sprintf( 'Empirical_omega_pdf_AM.fig' );
		saveas( fig16, fullfile( folderName, p_omega_emp_AM_filename ) );
		
		p_omega_emp_GD_filename = sprintf( 'Empirical_omega_pdf_GD.fig' );
		saveas( fig17, fullfile( folderName, p_omega_emp_GD_filename ) );
		
		pn_omega_optim_filename = sprintf( 'Optimal_omega_pdf.fig' );
		saveas( fig18, fullfile( folderName, pn_omega_optim_filename ) );
        
		V_2d_plot_status = 1;
		
	end
	
    function FT_integrand = v1_FT_func( x_1, x_2, omega_1, omega_2 )
		FT_integrand = v_1_func( x_1, x_2 ) .* exp( -1i * ( omega_1 * x_1 + omega_2 * x_2 ) ) / ( 2 * pi )^2;
	end
	
	
	function v1_value = v_1_func( x1, x2 )
		v1_value = ( x1.^2 / 2 + sqrt( 2 ) * x2.^2 / 2 + 2 * sin( x1 .* x2 ) .* chi_r_func( sqrt( x1.^2 + x2.^2 ), 2, 1 ) ) .* chi_r_func( sqrt( x1.^2 + x2.^2 ), 4, 0.2 );
    end