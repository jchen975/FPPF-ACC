function mpc = modify_case_load(casename, c)
	mpc = loadcase(casename);
	if nargin == 1
		c = 2.0;      % default behavior, scale up the load demand by 2.
	end
	mpc.bus(:, 3) = mpc.bus(:, 3) * c;
	mpc.bus(:, 4) = mpc.bus(:, 4) * c;
	
	fn = strcat('E:/Work/MASc/Research/Thesis/Nbus/Code/test_data/', casename, '_highload.m');
	fn_cycle = strcat(casename, '.mat');
	fn_cycle_new = strcat('E:/Work/MASc/Research/Thesis/Nbus/Code/cycle_matrix_data/', casename, '_highload_C.mat');
	if size(mpc.bus, 1) > 300 && exist(fn_cycle, 'file') == 2
		copyfile(fn_cycle, fn_cycle_new);
	end
	savecase(fn, mpc);
end