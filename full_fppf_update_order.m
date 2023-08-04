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
% case_name = 'case1888rte';   % voltage result doesn't match
% case_name = 'case1951rte';   
case_name = 'case2868rte';   % same here 
% case_name = 'case2869pegase';
% case_name = 'case6468rte'; 
% case_name = 'case9241pegase';

%% load test system
mpc = loadcase(case_name); 
mpc.case_name = case_name; 

%% modify r/x ratio as needed
max_rx = 0.8;
mpc = modify_case_rx_new(mpc, max_rx);

%% base/high load
load_type = 'high'; 
mpc.load_type = load_type;
mpc = modify_mpc_base_or_high_load(mpc, load_type);

%% flat start initial condition
mpc.bus(mpc.bus(:, 2)==1, 8) = 1;	% VM
mpc.bus(:, 9) = 0;					% VA

%% run full FPPF
% set FPPF options

% order = 'v-xc-psi';
% order = 'v-psi-xc'; 
% order = 'xc-psi-v';
% order = 'xc-v-psi';
% order = 'psi-v-xc';
% order = 'psi-xc-v'; 

tol = 1e-8;
max_iter = 1000;
orders = {'v-xc-psi', 'v-psi-xc', 'xc-psi-v', 'xc-v-psi', 'psi-v-xc', 'psi-xc-v'};

% run FPPF and extract results
for i = 1:length(orders)
	results = run_full_fppf_iterations(mpc, orders{i}, 'flat', tol, max_iter); 

	if results.success
		mismatch = results.mismatch;
		save(['test_data/mismatch/update_order/', case_name, '_order=', orders{i}, ...
			'_', load_type, '_load.mat'], 'mismatch');
	else
		fprintf('FPPF did not converge (%s, load type: %s, order: %s, ', case_name, load_type, orders{i})
		fprintf('psi-norm-big = %d). Did not save any data.\n', results.psi_norm_big);
	end
end

%% plot mismatches
plot_mismatch_update_order(case_name, load_type) 