function mpc = modify_mpc_base_or_high_load(mpc, load_type)
	if strcmp(load_type, 'base')
		gain = 1;  % base load
	elseif strcmp(load_type, 'high')
		gain = max(1, 0.9 * compute_lambda_2017a(mpc));  % high load - 90% of insolvability by NR
	end

	mpc.gen(:, 2) = gain*mpc.gen(:, 2);  % increased P generation
	mpc.bus(:, 3) = gain*mpc.bus(:, 3);  % and increased P load
	mpc.bus(:, 4) = gain*mpc.bus(:, 4);  % and increased Q load 
end