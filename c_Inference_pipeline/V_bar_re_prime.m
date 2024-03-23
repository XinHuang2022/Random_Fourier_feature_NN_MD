    
function V_bar_prime = V_bar_re_prime( eta_re, eta_im, omega, x, alpha, gamma, R_minus_2, R_minus_1, sigma_a, sigma_b )
		% x is a row vector of length 1 * dim
		% omega is a matrix of size dim * K
		
		% dim = size( omega, 1 );    % dimensionality is the number of rows in the omega array
		% K = length( eta_re );
		
		inner_prod = x * omega;
		dactivation_dx_part = -eta_re .* sin( inner_prod ) - eta_im .* cos( inner_prod );
		V_bar_trained_prime = dactivation_dx_part * omega';
		% dactivation_dx_double = [ dactivation_dx; dactivation_dx ];
		% gradient_sum = sum( dactivation_dx_double .* omega, 2 );

        % activation = eta_re .* cos( inner_prod ) - eta_im .* sin( inner_prod );
        % V_bar_re = sum( activation, 2 );
		V_bar_re = cos( inner_prod ) * eta_re' - sin( inner_prod ) * eta_im';
		
		% V_bar_trained_prime = gradient_sum';
		% r = norm( x );
		r = sqrt( x( :, 1 ).^2 + x( :, 2 ).^2 );
		% V_bar_prime_int = V_bar_trained_prime * chi_r_func( r, R_minus_1, R_0 ) + V_bar_re * dchi_dr_func( r, R_minus_1, R_0 ) / r * x;
		V_bar_prime_int = V_bar_trained_prime .* chi_r_func( r, R_minus_1, sigma_b ) + ( V_bar_re .* dchi_dr_func( r, R_minus_1, sigma_b ) ./ r ) .* x;
		
		V_value = V_func( x, alpha, gamma, R_minus_2, sigma_a );
		V_grad = V_prime( x, alpha, gamma, R_minus_2, sigma_a );
		chi_0_value = chi_r_func( r, R_minus_1, sigma_b );
		dchi_0_dr_value = dchi_dr_func( r, R_minus_1, sigma_b );
		V_bar_prime_ext = V_grad .* ( 1 - chi_0_value ) - ( V_value .* dchi_0_dr_value ./ r ) .* x;
		
		V_bar_prime = V_bar_prime_int + V_bar_prime_ext;
	
	end
	%{
	function gradient = compute_gradient( model, X )
		Y = forward( model, X );
		gradient = dlgradient( Y, X );
	end
	%}