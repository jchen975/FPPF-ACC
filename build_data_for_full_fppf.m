function data = build_data_for_full_fppf(mpc, verbose)

	case_name = mpc.case_name;
	if nargin == 1
		verbose = 0;
	end

	%% re-order and load system data
	
	[load_idx, gen_idx, gen_bus_idx_unique, edge_idx, ~, bus_idx_reverse, mpc_LG, A, Y, Yft, Ytf] = reorder_mpc_LG(mpc);   % mpc_LG is mpc with buses in the LG order
	idx_og2LG = [load_idx; gen_idx];	% mpc.bus(idx_og2LG, i) = mpc_LG.bus(:, i)
	idx_LG2og = bus_idx_reverse;		% mpc_LG.bus(idx_LG2og, i) = mpc.bus(:, i)
	% ALSO for conversion: mpc.bus(load_idx, 1) = mpc_LG.bus(1:n, 1), similar for gen

	% Get constants
	n = length(load_idx);
	m = length(gen_bus_idx_unique);
	N = n+m;    assert(N == size(mpc.bus,1));
	E = length(edge_idx);
	ref_idx = find(mpc_LG.bus(:,2) == 3);

	% Get real/reactive power injections; if there are multiple generators at a
	% bus, the gen bus PG is the sum of individual PG's.
	% [P, Q] = compute_power_injection(mpc_LG.bus, mpc_LG.gen, mpc_LG.baseMVA, n, m);
	% QL = Q(1:n);
	% !! TODO: check active generators and stop forcing all generators to be active in 
	% reorder_mpc_LG.m
	PD_load = mpc_LG.bus(1:n, 3);
	QD_load = mpc_LG.bus(1:n, 4);
	PD_gen  = mpc_LG.bus(n+1:end, 3);
	QD_gen  = mpc_LG.bus(n+1:end, 4);
	PG_gen  = mpc_LG.gen(gen_bus_idx_unique, 2);  %% FIX THIS!
	QG_gen  = mpc_LG.gen(gen_bus_idx_unique, 3);  %% AND THIS!

	P  = [-PD_load; (PG_gen - PD_gen)] ./ mpc_LG.baseMVA;
% 	Q  = [-QD_load; (QG_gen - QD_gen)] ./ mpc_LG.baseMVA;
	QL = -QD_load ./ mpc_LG.baseMVA;

	%% Open circuit voltage magnitude
	VG = mpc_LG.gen(gen_bus_idx_unique, 6);
	Vstar_L = compute_open_circuit_voltage(imag(Y), VG);
	Vstar = [Vstar_L; VG];

	%% Create/Compute matrices needed to run FPPF
	% Cycle matrix; if case has 300+ buses and previously computed C, then load it
	% otherwise compute again
	cycle_matrix_fn = strcat('cycle_matrix_data/', case_name, '_C.mat');
	if n+m >= 300 && isfile(cycle_matrix_fn) 
		if verbose
			fprintf('Found existing cycle matrix, loading time: ');
		end
		tic;
		C = load(cycle_matrix_fn, 'C').C;
		t = toc;
		if verbose
			fprintf('%.8f seconds.\n', t);
		end
	else
		if verbose
			fprintf('Did not find existing cycle matrix, computation time: ');
		end
		tic;
		C = sparse(null(full(A), 'r'));
		t = toc;
		if verbose
			fprintf('%.8f seconds.\n', t);
		end
		
		if n+m >= 300
			if verbose
				fprintf('Saving the computed cycle matrix, saving time: ');
			end
			tic;
			save(cycle_matrix_fn, 'C');
			t = toc;
			if verbose
				fprintf('%.8f seconds.\n', t);
			end
		end	
	end
% 	nc = size(C, 2);
	% assert(norm(full(A*C)) == 0);  % sanity check

	% Admittance (sub)matrix
	G = real(Y); B = imag(Y);
	Bdiag_LL = sparse(diag(diag(B(1:n,1:n))));
	Gdiag    = sparse(diag(diag(G)));
	BLL      = imag(Y(1:n, 1:n)); 

	% get the real and imaginary off-diagonal entries of Y, as arrays
	Gft = real(Yft);   Bft = imag(Yft);
	Gtf = real(Ytf);   Btf = imag(Ytf);

	% get from/to partition of the incidence matrix
	A_from = sparse(zeros(size(A))); 
	A_to   = sparse(zeros(size(A)));
	A_from(A == 1) = 1;
	A_to(A == -1) = 1;

	% the Edge x Edge [Vi*Vj*G], [Vi*Vj*B] matrices for both ij and ji 
	DG_ft = sparse(diag( Gft .* (A_from' * Vstar) .* (A_to' * Vstar) ));
	DG_tf = sparse(diag( Gtf .* (A_from' * Vstar) .* (A_to' * Vstar) ));
	DB_ft = sparse(diag( Bft .* (A_from' * Vstar) .* (A_to' * Vstar) ));
	DB_tf = sparse(diag( Btf .* (A_from' * Vstar) .* (A_to' * Vstar) )); 

	% the asymmetrically weighted incidence matrices, weighted by the diagonal
	% matrices defined above
	GammaG_directed   = A_from * DG_ft - A_to * DG_tf;
	GammaG_undirected = A_from * DG_ft + A_to * DG_tf;
	GammaB_directed   = A_from * DB_ft - A_to * DB_tf;
	GammaB_undirected = A_from * DB_ft + A_to * DB_tf;

	%%

	% stiffness matrix
	S = (1/4) * sparse(diag(Vstar_L)) * BLL * sparse(diag(Vstar_L));

	% R matrix for removing slack elements; FOR NOW only 1 slack bus @ ref
% 	alpha = zeros(n+m, 1);
% 	alpha(ref_idx) = 1.0;
	R = eye(n+m);
	R(:, ref_idx) = [];  % R is (n+m) x (n+m-1) and columns of R are orthonormal to alpha
	% assert(isequal(size(R), [n+m, n+m-1])); % sanity check
	% assert(norm(R'*alpha) == 0);			% sanity check
	Rt = sparse(R');

	% MB matrix (of size (n+m-1) x E), its pseudoinverse and its kernel matrix K

	MB = sparse(Rt * GammaB_directed);

	% use LUQ routine to compute K (method 2)
% 	fprintf('Using LUQ decomposition to find K matrix; computation time: ')
	[~, RNull] = spspaces(MB, 2);
	idx_luq = RNull{3};
	q_luq = RNull{1};
	K = q_luq(:, idx_luq);
	K = normalize(K, 'norm');  % normalize columns w/ 2-norm?? This has really large 1,2,inf norms

	MB_pinv = sparse(lsqminnorm(MB, eye(n+m-1)));

	data = struct('P', P, 'QL', QL, 'Vstar', Vstar, 'A', A, 'Rt', Rt, 'K', K, 'C', C, ...
		'S', S, 'MB', MB, 'MB_pinv', MB_pinv, 'Bdiag_LL', Bdiag_LL, 'Gdiag', Gdiag, ...
		 'GammaG_directed', GammaG_directed, 'GammaG_undirected', GammaG_undirected, ...
		 'GammaB_directed', GammaB_directed, 'GammaB_undirected', GammaB_undirected, ...
		 'idx_og2LG', idx_og2LG, 'idx_LG2og', idx_LG2og, 'mpc_LG', mpc_LG, 'ref_idx', ref_idx ...
	);
end