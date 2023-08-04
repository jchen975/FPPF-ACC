% clc
close all

current_dir = pwd;
if ~strcmp(current_dir, 'E:\Work\MASc\Research\Thesis\Nbus\Code')
	cd 'E:\Work\MASc\Research\Thesis\Nbus\Code'
	matlab_latex_setup;
end

clear

%% Test system to run
% case_name = 'C:/Program Files/MATLAB/R2021a/addons/matpower7.1/matpower7.1/data/case2';
% case_name = 'case2';
% case_name = 'case9';
case_name = 'case30';
% case_name = 'case89pegase';
% case_name = 'case118';


%% load test system
mpc = loadcase(case_name); 
mpc.case_name = case_name; 

%% base/high load
load_type = 'high'; 
mpc.load_type = load_type;
mpc = modify_mpc_base_or_high_load(mpc, load_type);

%% modify r/x ratio as needed
max_rx = 0.8;
mpc = modify_case_rx_new(mpc, max_rx);

%% run NR with good IC to get the known high voltage solution
mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.alg', 'NR', 'pf.tol', 1e-10, 'pf.nr.max_it', 100); 
ret_HV = runpf(mpc, mpopt);
if ~ret_HV.success
	input('This case is not solvable with the current IC. Press Ctrl+C to stop execution and try another IC.\n');
end
VM_HVsol = ret_HV.bus(:, 8);
VA_HVsol = ret_HV.bus(:, 9);

ref_idx = mpc.bus(mpc.bus(:, 2) == 3, 1);
if VA_HVsol(ref_idx) ~= 0	% ensure ref bus has 0 phase
	VA_HVsol = VA_HVsol - VA_HVsol(ref_idx);
end
%% run full FPPF
% set FPPF options

order = 'v-xc-psi';
% order = 'v-psi-xc'; 
% order = 'xc-psi-v';
% order = 'xc-v-psi';
% order = 'psi-v-xc';
% order = 'psi-xc-v'; 
tol = 1e-9;
max_iter = 1000;

%% sensitivity test
% verbose = 0;
% data = build_data_for_full_fppf(mpc, verbose);

fp_success = 0;		fp_big_psi = 0;
fd_success = 0;
nr_success = 0; 

spread = 0.1; 
ulim_v = 1+spread;		llim_v = 1-spread; 
ulim_t = 0;				llim_t = 0;
% ulim_v = 1;				llim_v = 1; 
% ulim_t = 45*spread;		llim_t = -45*spread;	


rng(99);
numIC = 1000;
tol2 = 1e-7;  % want a slightly more relaxed tolerance

% fprintf('Starting IC sensitivity test with %d samples.\n', numIC);
for i = 1:numIC
	% set IC
	n = length( mpc.bus(mpc.bus(:, 2) == 1, 1) );
	m = length( mpc.bus(:, 1) ) - n;
	
	VL_init = llim_v + (ulim_v-llim_v) .* rand(n, 1); 
    Va_init = llim_t + (ulim_t-llim_t) .* rand(n+m, 1);  
	mpc.bus(mpc.bus(:, 2) == 1, 8) = VL_init;
	mpc.bus(:, 9) = Va_init;
	
	% run FPPF
	results = run_full_fppf_iterations(mpc, order, '', tol, max_iter); 
	if results.success && norm(results.VMsol - VM_HVsol, Inf) <= tol2 && ...
						  norm(results.VAsol - VA_HVsol, Inf) <= tol2
		fp_success = fp_success + 1;
	elseif isfield(results, 'psi_norm_big')
		if results.psi_norm_big  % FPPF failed due to psi norm too big at some point
			fp_big_psi = fp_big_psi + 1;
		end
	end
	
	% Compare against NR
	mpopt_NR = mpoption('verbose', 0, 'out.all', 0, 'pf.alg', 'NR', 'pf.tol', tol, 'pf.nr.max_it', max_iter);
	ret_nr = runpf(mpc, mpopt_NR); 
	if ret_nr.bus(ref_idx,9) ~= 0	% ensure ref bus has 0 phase
		ret_nr.bus(:,9) = ret_nr.bus(:,9) - ret_nr.bus(ref_idx,9);
	end
	if ret_nr.success && norm(ret_nr.bus(:,8) - VM_HVsol, Inf) <= tol2 && ...
						 norm(ret_nr.bus(:,9) - VA_HVsol, Inf) <= tol2
		nr_success = nr_success + 1; 
	end
	
	% Compare against FDLF
	mpopt_NR = mpoption('verbose', 0, 'out.all', 0, 'pf.alg', 'FDXB', 'pf.tol', tol, 'pf.fd.max_it', max_iter);
	ret_fd = runpf(mpc, mpopt_NR); 
	if ret_fd.bus(ref_idx,9) ~= 0	% ensure ref bus has 0 phase
		ret_fd.bus(:,9) = ret_fd.bus(:,9) - ret_fd.bus(ref_idx,9);
	end
	if ret_fd.success && norm(ret_fd.bus(:,8) - VM_HVsol, Inf) <= tol2 && ...
						 norm(ret_fd.bus(:,9) - VA_HVsol, Inf) <= tol2
		fd_success = fd_success + 1; 
	end
end

fprintf('Initialization sensitivity test (%d samples, spread=%.2f, load type=%s, order=%s) complete: \n', numIC, spread, mpc.load_type, order); 
fprintf(' - NR   converged %d times (%.2f %%)\n', nr_success, nr_success/numIC * 100); 
fprintf(' - FDLF converged %d times (%.2f %%)\n', fd_success, fd_success/numIC * 100);
fprintf(' - FPPF converged %d times (%.2f %%)\n', fp_success, fp_success/numIC * 100);
if fp_success < numIC
	fprintf('   -> FPPF did not converge in %d cases (%.2f %%), %d of them due to big psi norm\n', ...
		numIC-fp_success, (numIC-fp_success)/numIC * 100, fp_big_psi);
end
