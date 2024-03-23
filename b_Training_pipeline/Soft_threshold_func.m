    
	function S_gamma_value = Soft_threshold_func( u, gamma, K )
		
		S_gamma_value = ( u > gamma ) .* ( u - gamma ) + ( u < -gamma ) .* ( u + gamma );
        
	end
		
		
		
		