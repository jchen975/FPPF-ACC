function mpc = modify_case_rx(casename, ulim, N)
	rng(0);
	mpc = loadcase(casename);
	if nargin == 1
		ulim = 0.9;		% default upper limit on branch R/X ratio is 0.9
	elseif nargin == 2
		N = 0;          % if 0, cap the high R/X ratio to be ulim; otherwise, change N random lines R/X ratios 
	end
	branch = mpc.branch;
	
	found_large_rx = false;
	if N == 0
		for i = 1:size(branch(:,1))
			rx_i = branch(i, 3) / branch(i, 4);
			if rx_i > ulim
				found_large_rx = true; 
				prompt = strcat('Branch (', num2str(branch(i,1)), ', ', num2str(branch(i,2)), ...
						') has R/X ratio: ', num2str(rx_i), ', greater than the specified upper limit: ', ...
						num2str(ulim), '. Type 1 to change to the upper limit, and 2 to ignore.');
				answer = input(prompt); 
				if answer ~= 1
					continue;
				elseif answer == 1
					branch(i, 3) = ulim * branch(i, 4);
					rx_i_new = branch(i, 3) / branch(i, 4);
					fprintf('>>> R/X ratio of branch (%4d, %4d) changed from %5.4f to %5.4f.\n', branch(i,1), branch(i,2), rx_i, rx_i_new);
				end
			end
		end
		if ~found_large_rx
			fprintf('No branch R/X ratio exceeding %.3f exists in this system.\n', ulim);
		end
	else 
		E = size(mpc.branch, 1);
		idx = randi(E, N, 1);    % change N random branches to have high RX ratio = ulim
		for i = 1:N
			rx_i = branch(idx(i), 3) / branch(idx(i), 4);
			mpc.branch(idx(i), 3) = mpc.branch(idx(i), 4) * ulim;
			fprintf('R/X ratio of branch (%4d, %4d) changed from %5.4f to %5.4f.\n', branch(idx(i),1), branch(idx(i),2), rx_i, ulim);
		end
	end
	
	mpc.branch = branch;
	
	if ulim <= 1.0
		ftype = '_lowRX';
	else
		ftype = '_highRX';
	end
	
	fn = strcat('E:/Work/MASc/Research/Thesis/Nbus/Code/test_data/', casename, ftype, '.m');
	fn_cycle = strcat(casename, '.mat');
	if size(mpc.bus, 1) > 300 && exist(fn_cycle, 'file') == 2
		fn_cycle_new = strcat('E:/Work/MASc/Research/Thesis/Nbus/Code/cycle_matrix_data/', casename, ftype, '_C.mat');
		copyfile(fn_cycle, fn_cycle_new);
	end
	if found_large_rx || N > 0   % save the modified case file
		savecase(fn, mpc);
	end
end