    
	
Num_replica = 32;
M_sp = 2^16;

tau = 2;
Delta_tau_plot = 0.1;

Num_parallel_worker = 8;
L_x_vectorize = 16;
h_x_sample = 0.005;

filename0 = 'testrun_true_V_correlation';
datadir0 = './';


exitcode = test_approx_V_md_corr( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, Num_parallel_worker )
	
	
	