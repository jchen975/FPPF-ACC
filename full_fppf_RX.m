clc
close all

current_dir = pwd;
if ~strcmp(current_dir, 'E:\Work\MASc\Research\Thesis\Nbus\Code')
	cd 'E:\Work\MASc\Research\Thesis\Nbus\Code'
	matlab_latex_setup;
end
% clear

%% test system to run

% case_name = 'C:/Program Files/MATLAB/R2021a/addons/matpower7.1/matpower7.1/data/case2';
% case_name = 'case2';
% case_name = 'case9';
% case_name = 'case30';
% case_name = 'case89pegase';
% case_name = 'case118';
% case_name = 'case300';
% case_name = 'case1354pegase';
case_name = 'case2869pegase';
% case_name = 'case6468rte'; 
% case_name = 'case9241pegase';

%% load test system
mpc = loadcase(case_name); 
mpc.case_name = case_name; 

%% base/high load
load_type = 'base'; 
mpc.load_type = load_type;
mpc = modify_mpc_base_or_high_load(mpc, load_type);

%% flat start initial condition
mpc.bus(mpc.bus(:, 2)==1, 8) = 1;	% VM
mpc.bus(:, 9) = 0;					% VA

%% run full FPPF
% set FPPF options
order = 'v-xc-psi';
% order = 'xc-psi-v';
% order = 'xc-v-psi';
% order = 'v-psi-xc'; 
% order = 'psi-v-xc';
% order = 'psi-xc-v';

tol = 1e-8;
max_iter = 1000;

%% run FPPF with different r/x ratios 
max_rx = [0.5, 0.6 0.7, 0.8, 0.85, 0.9, 0.95, 0.99];  % case300

for i = 1:length(max_rx)
	mpc_copy = mpc;
	mpc_copy = modify_case_rx_new(mpc_copy, max_rx(i));

	% run FPPF and extract results
	results = run_full_fppf_iterations(mpc_copy, order, 'flat', tol, max_iter); 

	if results.success
		mismatch = results.mismatch;
		save(['test_data/mismatch/rx/', case_name, '_rx=', num2str(max_rx(i),2), ...
				'_', load_type, '_load_order=', order, '.mat'], 'mismatch');
	else
		fprintf('FPPF did not converge (order: %s, psi-norm-big = %d). Did not save any data.\n', order, results.psi_norm_big);
	end
end

%%
plot_mismatch_RX(case_name, max_rx, load_type);