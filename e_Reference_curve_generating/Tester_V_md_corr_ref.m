%
% Tester for True_V_md_corr_run.m function
%

save_ref_curve = 1;
	
Num_replica = 256;
M_sp = 2^20;

tau = 6;
Delta_tau_plot = 0.1;

Num_parallel_worker = 128;
L_x_vectorize = 16;
h_x_sample = 0.005;

filename0 = 'testrun_true_V_correlation';
datadir0 = './';

%{
filename0='testrun_true_V_correlation';
lulog0='1';
datadir0='./';
dim0 = '3';
Ne0 = '6';
beta0 = '1';
lambda0 = '0';
dt0 = '0.025';
Mxmax0 = '18';
dMx0 = '0.5';
Nrep0 = '10';
Ncore0 = '128';
TestCase0 = "New_V1";
%}

exitcode = True_V_md_corr_run( h_x_sample, M_sp, Num_replica, tau, Delta_tau_plot, ...
							   Num_parallel_worker, L_x_vectorize, save_ref_curve, filename0, datadir0 )

%%%
% dim0:      the dimension of the space.
% Ne0:       the number of electrons.
% beta0:     the inverse temperature 1 / ( k_B * T ).
% lambda0:   the relative force parameter of Coulomb interactions. 
%            For the fermi gas case under Harmonic Oscillator potential 
%            without Coulomb interactions, lambda will be automatically 
%            set to 0.
% dt0:       the time step for the integration in the exponent of the 
%            Feymann-Kac representation.
% 2^Mxmax0:  the sample size for each independent estimator of partition 
%            function Z and mean-field h.
% dMx0:      the x-axis distance between two points in the figures for the 
%            convergence of estimators with increasing sample size.
% 2^Nrep0:   the total number of independent replicas of 
%            partition function Z and mean-field h, each replica generates
%            and stores 2^Mxmax0 samples for Monte Carlo integral.
% Ncore0:    the number of parallel workers for evaluating all the 
%            2^Nrep0 replicas.
% TestCase0: Case V1 denotes Harmonic Oscillator potential, Case V2 for 
%            Coulomb interaction under a harmonic trap.
