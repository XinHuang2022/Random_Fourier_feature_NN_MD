     
    clear all
    close all
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
    
	alpha = sqrt( 2 );
	gamma = 2;
	
    R_a = 2;
	R_b = 3;
	R_c = 4;
    sigma_a = 1;
	sigma_b = 0.2;
	sigma_c = 0.2;
	
	% Compute the function values for the grid of (x, y) points
	R = sqrt( X.^2 + Y.^2 );
	Z_true = V_true_2d_func( X, Y, alpha, gamma, R_a, sigma_a ) .* chi_r_func( R, R_b, sigma_b ) ;
	% Z_approx = V_approx_2d_plot( X, Y, R_0, R_1, alpha_chi );

	% Create a 3D surface plot
	figure( 1 );
	surf( X, Y, Z_true, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );

	% Customize the plot
	xlabel( '$x_1$-axis', 'interpreter', 'latex' );
	ylabel( '$x_2$-axis', 'interpreter', 'latex' );
	zlabel( '$V_i(x_1, x_2)$', 'interpreter', 'latex' );
	title( '3D Plot of the Function $V_i = \big(\frac{ x_1^2}{2} + \frac{ \alpha x_2^2}{2} + \gamma \sin{( x_1 x_2 )} * \chi_{a}(r)\big)\chi_b(r)$', 'interpreter', 'latex');

	% You can also customize the colormap and add a color bar if needed
	colormap('jet'); % Change the colormap
	colorbar;        % Add a color bar

	% Adjust the view angle if necessary
	view(0, 90); % Specify the view angle (azimuth, elevation)
	
	
	Vi_grad_x1 = Vi_grad_x1_func( X, Y, alpha, gamma, R_a, sigma_a, R_b, sigma_b );
	Vi_grad_x2 = Vi_grad_x2_func( X, Y, alpha, gamma, R_a, sigma_a, R_b, sigma_b );

	figure( 2 );
	surf( X, Y, Vi_grad_x1, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );
	xlabel( '$x_1$-axis', 'interpreter', 'latex' );
	ylabel( '$x_2$-axis', 'interpreter', 'latex' );
	zlabel( '$\nabla_{x_1} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	title( 'Plot of the derivative of the function $\nabla_{x_1} V(x_1, x_2)$', 'interpreter', 'latex' );

	colormap('jet'); % Change the colormap
	colorbar;        % Add a color bar

	view(0, 90); % Specify the view angle (azimuth, elevation)
	
	
	figure( 3 );
	surf( X, Y, Vi_grad_x2, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );
	xlabel( '$x_1$-axis', 'interpreter', 'latex' );
	ylabel( '$x_2$-axis', 'interpreter', 'latex' );
	zlabel( '$\nabla_{x_2} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	title( 'Plot of the derivative of the function $\nabla_{x_2} V(x_1, x_2)$', 'interpreter', 'latex' );

	colormap('jet'); % Change the colormap
	colorbar;        % Add a color bar

	view(0, 90); % Specify the view angle (azimuth, elevation)

    
    Vi_Hessian_11 = Vi_Hess_11_func( X, Y, alpha, gamma, R_a, sigma_a, R_b, sigma_b );
	figure( 4 );
	surf( X, Y, Vi_Hessian_11, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );
	xlabel( '$x_1$-axis', 'interpreter', 'latex' );
	ylabel( '$x_2$-axis', 'interpreter', 'latex' );
	zlabel( '$\frac{\partial^2}{\partial{x_1}^2} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	title( 'Plot of the derivative of the function $\frac{\partial^2}{\partial{x_1}^2} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	
	Vi_Hessian_12 = Vi_Hess_12_func( X, Y, alpha, gamma, R_a, sigma_a, R_b, sigma_b );
	figure( 5 );
	surf( X, Y, Vi_Hessian_12, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );
	xlabel( '$x_1$-axis', 'interpreter', 'latex' );
	ylabel( '$x_2$-axis', 'interpreter', 'latex' );
	zlabel( '$\frac{\partial^2}{\partial{x_1}\partial{x_2}} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	title( 'Plot of the derivative of the function $\frac{\partial^2}{\partial{x_1}\partial{x_2}} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	
	Vi_Hessian_22 = Vi_Hess_22_func( X, Y, alpha, gamma, R_a, sigma_a, R_b, sigma_b );
	figure( 6 );
	surf( X, Y, Vi_Hessian_22, 'EdgeColor', 'none' );
    % contour( X, Y, Z_true, 'EdgeColor', 'none' );
	xlabel( '$x_1$-axis', 'interpreter', 'latex' );
	ylabel( '$x_2$-axis', 'interpreter', 'latex' );
	zlabel( '$\frac{\partial^2}{\partial_{x_2}^2} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	title( 'Plot of the derivative of the function $\frac{\partial^2}{\partial{x_2}^2} V_i(x_1, x_2)$', 'interpreter', 'latex' );
	
	
    
	%{
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
	%}

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
	