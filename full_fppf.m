% clc
% close all

current_dir = pwd;
if ~strcmp(current_dir, 'E:\Work\MASc\Research\Thesis\Nbus\Code')
	cd 'E:\Work\MASc\Research\Thesis\Nbus\Code'
	matlab_latex_setup;
end
clear

%% test system to run

% case_name = 'C:/Program Files/MATLAB/R2021a/addons/matpower7.1/matpower7.1/data/case2';
% case_name = 'case9';
% case_name = 'case30';  % 
% case_name = 'case89pegase';
case_name = 'case118';
% case_name = 'case300';
% case_name = 'case1354pegase';
% case_name = 'case1888rte';   % voltage result doesn't match
% case_name = 'case1951rte';   
% case_name = 'case2868rte';  % same here
% case_name = 'case2869pegase';
% case_name = 'case6468rte'; 
% case_name = 'case9241pegase';

%% load test system
mpc = loadcase(case_name); 
mpc.case_name = case_name; 

%% base/high load
load_type = 'base'; 
mpc.load_type = load_type;
mpc = modify_mpc_base_or_high_load(mpc, load_type);

%% modify r/x ratio as needed
max_rx = 0.8;
mpc = modify_case_rx_new(mpc, max_rx);

%% flat start initial condition
mpc.bus(mpc.bus(:, 2)==1, 8) = 1;	% VM
mpc.bus(:, 9) = 0;					% VA

%% run full FPPF
% set FPPF options
order = 'v-xc-psi';
% order = 'xc-psi-v';
% order = 'xc-v-psi';

% order='v-psi-xc'; 
% order='psi-v-xc';
% order='psi-xc-v';

tol = 1e-8;
max_iter = 100;

% run FPPF and extract results
results = run_full_fppf_iterations(mpc, order, 'flat', tol, max_iter); 

if results.success
	mismatch = results.mismatch;
% 	save(['test_data/mismatch/iterations/', case_name,'_', load_type, ... 
% 				'_load_order=', order, '.mat'], 'mismatch');
else
	fprintf('FPPF did not converge (%s load, order: %s, psi-norm-big = %d). Did not save any data.\n', load_type, order, results.psi_norm_big);
end


%% compare against NR

mpopt_NR = mpoption('verbose', 0, 'out.all', 0, 'pf.alg', 'NR', 'pf.tol', tol, 'pf.nr.max_it', max_iter);
ret_NR = runpf(mpc, mpopt_NR);
if ret_NR.success
end

mpopt_FD = mpoption('verbose', 0, 'out.all', 0, 'pf.alg', 'FDXB', 'pf.tol', tol, 'pf.fd.max_it', max_iter);
ret_FD = runpf(mpc, mpopt_FD);
assert(ret_NR.success == 1);

% N = length(results.VMsol);

% figure
% plot(1:N, results.VMsol, 'o');
% hold on
% plot(1:N, ret_NR.bus(:,8), '*');
% plot(1:N, ret_FD.bus(:,8), 's');
% legend('FPPF', 'NR', 'FDLF')
% title('VM');
% 
% figure
% plot(1:N, results.VAsol, 'o');
% hold on
% plot(1:N, ret_NR.bus(:,9), '*')
% plot(1:N, ret_FD.bus(:,9), 's')
% legend('FPPF', 'NR', 'FDLF')
% title('VA');