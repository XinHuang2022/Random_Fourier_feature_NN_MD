
	% Define the function you want to plot
	% In this example, we'll use a simple quadratic function as an example.
	% You should replace this with your actual function.
	% z = f(x, y)
	% f = @(x, y) x.^2 + y.^2;

	% Create a grid of x and y values
	% x = linspace(-2, 2, 100); % Define a range of x values
	% y = linspace(-2, 2, 100); % Define a range of y values
	x = -5 : 0.02 : 5;
	y = -5 : 0.02 : 5;
	[ X, Y ] = meshgrid( x, y );  % Create a grid of (x, y) values
    
    R_a = 2;
	R_b = 3;
	
	sigma_a = 1;
	sigma_b = 0.2;
	
	% Compute the function values for the grid of (x, y) points
	Z_true = V_true_2d_plot( X, Y, R_a, sigma_a );
	Z_approx = V_approx_2d_plot( X, Y, R_a, R_b, sigma_a, sigma_b );

	% Create a 3D surface plot
	figure( 1 );
	surf( X, Y, Z_true, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );

	% Customize the plot
	xlabel( 'x_1-axis' );
	ylabel( 'x_2-axis' );
	zlabel( 'V(x_1, x_2)' );
	title( '3D Plot of the Function V = 0.5 x_1^2 + 0.5 \alpha x_2^2 + \gamma sin( x_1 x_2 ) * \chi_{-1}(r)');

	% You can also customize the colormap and add a color bar if needed
	colormap('jet'); % Change the colormap
	colorbar;        % Add a color bar

	% Adjust the view angle if necessary
	view(0, 90); % Specify the view angle (azimuth, elevation)
    

    figure( 2 );
	surf( X, Y, Z_approx, 'EdgeColor', 'none' );
    % Customize the plot
	xlabel( 'x_1-axis' );
	ylabel( 'x_2-axis' );
	zlabel( '\bar{V}(x_1, x_2)' );
	title( '3D Plot of the approximated function \bar{V}');

	% You can also customize the colormap and add a color bar if needed
	colormap('jet'); % Change the colormap
	colorbar;        % Add a color bar

	% Adjust the view angle if necessary
	view(0, 90); % Specify the view angle (azimuth, elevation)


    figure( 3 );
	surf( X, Y, abs( Z_approx - Z_true ), 'EdgeColor', 'none' );
    % Customize the plot
	xlabel( 'x_1-axis' );
	ylabel( 'x_2-axis' );
	zlabel( '|V(x_1, x_2) - \bar{V}(x_1, x_2)|' );
	title( '3D Plot of the difference between V(x) and \bar{V}(x)');

	% You can also customize the colormap and add a color bar if needed
	colormap('jet'); % Change the colormap
	colorbar;        % Add a color bar

	% Adjust the view angle if necessary
	view(0, 90); % Specify the view angle (azimuth, elevation)


	%{
	% Adjust the axis limits for a better view
	xlim([-2, 2]);
	ylim([-2, 2]);
	zlim([0, 8]);
    %}
    %{
    figure;
    plot( x, Z( 251, : ) );
    xlabel( 'x_1-axis' );
	ylabel( 'V(x)' );
    title( 'Section of the potential V(x_1,x_2) on x_2 = 1' )
    %}
	