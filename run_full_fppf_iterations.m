function result = run_full_fppf_iterations(data, order, IC, tol, max_iter, verbose) 
% 	INPUTS: 
% 		- data: matpower system struct (mpc) or data from
%				build_data_for_full_fppf.m (if modify R/X ratio, do it before here!)
% 		- order: string for update order of the variables separated by `-`, 
% 			e.g., 'v-xc-psi' updates v, then xc, then psi
% 		- IC: initial condition type, if 'flat' then run with flat start;
%			otherwise run with default mpc values
%		- tol: power balance mistmatch tolerance
% 		- max_iter: maximum iterations for FPPF
%		- verbose: output verbosity
% 		
% 	OUTPUT
% 		- result: a struct with the following fields
% 			- success: flag indicating whether FPPF converges
% 			- v: final v values
% 			- psi: final psi values
%			- VMsol: solution for Vm in original indexing
% 			- VAsol: solution for Va in original indexing
% 			- mismatch: mismatch vector
% 			- iterations: number of iterations
%		   (- psi_norm_big: 1 if true, 0 if did not converge for other reasons.
%		   only available if success == 0)

	if nargin == 5
		verbose = 0;
	end
	
	%% unpack data
	if isfield(data, 'bus')	 % data is mpc and need to compute the actual data
		data = build_data_for_full_fppf(data, verbose);
	end

	[P, QL, Vstar, A, Rt, K, C, S, MB_pinv, Bdiag_LL, Gdiag] = ... 
		deal(data.P, data.QL, data.Vstar, data.A, data.Rt, data.K, data.C, data.S, ...
	data.MB_pinv, data.Bdiag_LL, data.Gdiag);

	[GammaG_directed, GammaG_undirected, GammaB_directed, GammaB_undirected] = ...
		deal( data.GammaG_directed, ... 
	data.GammaG_undirected, data.GammaB_directed, data.GammaB_undirected);

	[idx_og2LG, idx_LG2og, mpc_LG, ref_idx] = deal(data.idx_og2LG, data.idx_LG2og, ...
		data.mpc_LG, data.ref_idx); 

	n = length(QL);
	m = length(P) - n;
	N = n + m;
	E = size(C, 1);
	nc = size(C, 2);

	Vstar_L = Vstar(1:n);
	VG = Vstar(n+1:n+m);
	load_idx_og2LG = idx_og2LG(1:n);
	gen_idx_og2LG  = idx_og2LG(n+1:n+m);
	
	%% Initialize Full FPPF

	psi      = zeros(E, max_iter);
	v        = zeros(n, max_iter);
	xc       = zeros(nc, max_iter);
	mismatch = zeros(1, max_iter);

	i = 1;
	if strcmp(IC, 'flat')		% flat start
		v(:,i)   = 1 ./ Vstar_L;   
	else						% default values
		psi(:,i) = sin(A'*deg2rad(mpc_LG.bus(:, 9)));
		v(:,i)   = mpc_LG.bus(1:n, 8) ./ Vstar_L; 
	end
	
	% Compute initial mismatch
	hv = h(v(:,1), A);
	mismatch(i) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
								GammaG_undirected, GammaB_directed, ...
								GammaB_undirected, v(:,1), psi(:,1), Vstar, A, Rt, hv);


	%% Run Full FPPF
	psi_norm_big = false;
	if verbose
		fprintf('\nRunning Full FPPF (%s load, %s update order) now.\n', mpc.load_type, order); 
	end
	if strcmp(order, 'xc-psi-v')
		while mismatch(i) > tol && i <= max_iter
			tic;
			gv = g(v(:,i), m);						% this is a (n,1) vector   
			hv = h(v(:,i), A);						% this is a (E,1) vector
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % same here
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!

			% update c
			if nc > 0
				J_cycle = sparse(C' * diag(ones(E,1) ./ sqrt_term_psi) * hv_inv * K);
				xc(:, i+1) = xc(:, i) - J_cycle \ (C' * asin(psi(:, i))); 
			end 

			% update psi
			p1 = Gdiag * (Vstar .* gv) .^2;
			p2 = GammaG_undirected * (hv .* sqrt_term_psi);
			soln_particular  = hv_inv * MB_pinv * Rt * (P - p1 - p2);
			if nc > 0
				soln_homogeneous = hv_inv * K * xc(:, i+1);        % issue here if network is radial
			else
				soln_homogeneous = zeros(size(soln_particular));
			end

			psi(:, i+1) = soln_particular + soln_homogeneous;

			if norm(psi(:, i+1), inf) > 1
				if verbose
					fprintf('Warning: inf-norm of psi is too large. Exiting the loop.\n');
				end
				psi_norm_big = true;
				break;
			end

			% update v 
			sqrt_term_psi = ones(E,1) - sqrt(ones(E,1) - psi(:,i+1).^2);
			q1 = GammaG_directed(1:n, :) * (hv .* psi(:,i+1));
			q2 = GammaB_undirected(1:n, :) * (hv .* sqrt_term_psi);
			v(:, i+1) = ones(n,1) - ((4*S) \ ((QL - q1 - q2) ./ v(:,i))); 

			% compute iteration mismatch
			mismatch(i+1) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
									GammaG_undirected, GammaB_directed, ...
									GammaB_undirected, v(:,i+1), psi(:,i+1), Vstar, A, Rt, hv);
			i = i+1;
		end
	elseif strcmp(order, 'v-xc-psi')
		while mismatch(i) > tol && i <= max_iter 
			hv = h(v(:,i), A);						% this is a (E,1) vector
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % same here

			% update v 
			q1 = GammaG_directed(1:n, :) * (hv .* psi(:,i));
			q2 = GammaB_undirected(1:n, :) * (hv .* (1-sqrt_term_psi));
			v(:, i+1) = ones(n,1) - ((4*S) \ ((QL - q1 - q2) ./ v(:,i)));

			hv = h(v(:,i+1), A);					% this is a (E,1) vector
			gv = g(v(:,i+1), m);					% this is a (n,1) vector   
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!

			% update xc
			if nc > 0
				J_cycle = C' * sparse(diag(ones(E,1) ./ sqrt_term_psi)) * hv_inv * K;
				xc(:, i+1) = xc(:, i) - J_cycle \ (C' * asin(psi(:, i))); 
			end

			% update psi
			p1 = Gdiag * (Vstar .* gv) .^2;
			p2 = GammaG_undirected * (hv .* sqrt_term_psi);
			soln_particular  = hv_inv * MB_pinv * Rt * (P - p1 - p2);
			if nc > 0
				soln_homogeneous = hv_inv * K * xc(:, i+1);
			else
				soln_homogeneous = 0; %zeros(size(soln_particular));
			end

			psi(:, i+1) = soln_particular + soln_homogeneous;

			if norm(psi(:, i+1), inf) > 1
				if verbose
					fprintf('Warning: inf-norm of psi is too large. Exiting the loop.\n');
				end
				psi_norm_big = true;
				break;
			end

			% compute iteration mismatch
			mismatch(i+1) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
									GammaG_undirected, GammaB_directed, ...
									GammaB_undirected, v(:,i+1), psi(:,i+1), Vstar, A, Rt, hv);
			i = i+1;
		end
	elseif strcmp(order, 'xc-v-psi')
		while mismatch(i) > tol && i <= max_iter 
			hv = h(v(:,i), A);						% this is a (E,1) vector
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % same here
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!
			
			% update xc
			if nc > 0
				J_cycle = C' * sparse(diag(ones(E,1) ./ sqrt_term_psi)) * hv_inv * K;
				xc(:, i+1) = xc(:, i) - J_cycle \ (C' * asin(psi(:, i))); 
			end
			
			% update v 
			q1 = GammaG_directed(1:n, :) * (hv .* psi(:,i));
			q2 = GammaB_undirected(1:n, :) * (hv .* (1-sqrt_term_psi));
			v(:, i+1) = ones(n,1) - ((4*S) \ ((QL - q1 - q2) ./ v(:,i)));

			hv = h(v(:,i+1), A);					% this is a (E,1) vector
			gv = g(v(:,i+1), m);					% this is a (n,1) vector   
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % this is a (E,1) vector

			% update psi
			p1 = Gdiag * (Vstar .* gv) .^2;
			p2 = GammaG_undirected * (hv .* sqrt_term_psi);
			soln_particular  = hv_inv * MB_pinv * Rt * (P - p1 - p2);
			if nc > 0
				soln_homogeneous = hv_inv * K * xc(:, i+1);
			else
				soln_homogeneous = 0; %zeros(size(soln_particular));
			end

			psi(:, i+1) = soln_particular + soln_homogeneous;

			if norm(psi(:, i+1), inf) > 1
				if verbose
					fprintf('Warning: inf-norm of psi is too large. Exiting the loop.\n');
				end
				psi_norm_big = true;
				break;
			end

			% compute iteration mismatch
			mismatch(i+1) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
									GammaG_undirected, GammaB_directed, ...
									GammaB_undirected, v(:,i+1), psi(:,i+1), Vstar, A, Rt, hv);
			i = i+1;
		end
	elseif strcmp(order, 'v-psi-xc')
		while mismatch(i) > tol && i <= max_iter 
			hv = h(v(:,i), A);						% this is a (E,1) vector
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % same here

			% update v 
			q1 = GammaG_directed(1:n, :) * (hv .* psi(:,i));
			q2 = GammaB_undirected(1:n, :) * (hv .* (1-sqrt_term_psi));
			v(:, i+1) = ones(n,1) - ((4*S) \ ((QL - q1 - q2) ./ v(:,i)));

			hv = h(v(:,i+1), A);					% this is a (E,1) vector
			gv = g(v(:,i+1), m);					% this is a (n,1) vector   
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!
			
			% update psi
			p1 = Gdiag * (Vstar .* gv) .^2;
			p2 = GammaG_undirected * (hv .* sqrt_term_psi);
			soln_particular  = hv_inv * MB_pinv * Rt * (P - p1 - p2);
			if nc > 0
				soln_homogeneous = hv_inv * K * xc(:, i);
			else
				soln_homogeneous = 0; %zeros(size(soln_particular));
			end

			psi(:, i+1) = soln_particular + soln_homogeneous;

			if norm(psi(:, i+1), inf) > 1
				if verbose
					fprintf('Warning: inf-norm of psi is too large. Exiting the loop.\n');
				end
				psi_norm_big = true;
				break;
			end

			% update xc
			if nc > 0
				J_cycle = C' * sparse(diag(ones(E,1) ./ sqrt_term_psi)) * hv_inv * K;
				xc(:, i+1) = xc(:, i) - J_cycle \ (C' * asin(psi(:, i+1))); 
			end

			% compute iteration mismatch
			mismatch(i+1) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
									GammaG_undirected, GammaB_directed, ...
									GammaB_undirected, v(:,i+1), psi(:,i+1), Vstar, A, Rt, hv);
			i = i+1;
		end
	elseif strcmp(order, 'psi-v-xc')
		while mismatch(i) > tol && i <= max_iter 
			hv = h(v(:,i), A);						% this is a (E,1) vector
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % same here
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!
			gv = g(v(:,i), m);					% this is a (n,1) vector   

			% update psi
			p1 = Gdiag * (Vstar .* gv) .^2;
			p2 = GammaG_undirected * (hv .* sqrt_term_psi);
			soln_particular  = hv_inv * MB_pinv * Rt * (P - p1 - p2);
			if nc > 0
				soln_homogeneous = hv_inv * K * xc(:, i);
			else
				soln_homogeneous = 0; %zeros(size(soln_particular));
			end

			psi(:, i+1) = soln_particular + soln_homogeneous;

			if norm(psi(:, i+1), inf) > 1
				if verbose
					fprintf('Warning: inf-norm of psi is too large. Exiting the loop.\n');
				end
				psi_norm_big = true;
				break;
			end
			
			sqrt_term_psi = sqrt(1 - psi(:,i+1).^2); 
			% update v
			q1 = GammaG_directed(1:n, :) * (hv .* psi(:,i+1));
			q2 = GammaB_undirected(1:n, :) * (hv .* (1-sqrt_term_psi));
			v(:, i+1) = ones(n,1) - ((4*S) \ ((QL - q1 - q2) ./ v(:,i)));

			hv = h(v(:,i+1), A);					% this is a (E,1) vector
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!

			% update xc
			if nc > 0
				J_cycle = C' * sparse(diag(ones(E,1) ./ sqrt_term_psi)) * hv_inv * K;
				xc(:, i+1) = xc(:, i) - J_cycle \ (C' * asin(psi(:, i+1))); 
			end

			% compute iteration mismatch
			mismatch(i+1) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
									GammaG_undirected, GammaB_directed, ...
									GammaB_undirected, v(:,i+1), psi(:,i+1), Vstar, A, Rt, hv);
			i = i+1;
		end
	elseif strcmp(order,'psi-xc-v')
		while mismatch(i) > tol && i <= max_iter 
			hv = h(v(:,i), A);						% this is a (E,1) vector
			sqrt_term_psi = sqrt(1 - psi(:,i).^2);  % same here
			gv = g(v(:,i), m);					    % this is a (n,1) vector   
			hv_inv = sparse(diag(1 ./ hv));         % NEED TO SPARSIFY HERE OR ELSE THE PARTICULAR SOLUTION TAKES FOREVER!
			
			% update psi
			p1 = Gdiag * (Vstar .* gv) .^2;
			p2 = GammaG_undirected * (hv .* sqrt_term_psi);
			soln_particular  = hv_inv * MB_pinv * Rt * (P - p1 - p2);
			if nc > 0
				soln_homogeneous = hv_inv * K * xc(:, i);
			else
				soln_homogeneous = 0; %zeros(size(soln_particular));
			end

			psi(:, i+1) = soln_particular + soln_homogeneous;

			if norm(psi(:, i+1), inf) > 1
				if verbose
					fprintf('Warning: inf-norm of psi is too large. Exiting the loop.\n');
				end
				psi_norm_big = true;
				break;
			end

			% update xc
			sqrt_term_psi = sqrt(1 - psi(:,i+1).^2);  % same here
			if nc > 0
				J_cycle = C' * sparse(diag(ones(E,1) ./ sqrt_term_psi)) * hv_inv * K;
				xc(:, i+1) = xc(:, i) - J_cycle \ (C' * asin(psi(:, i+1))); 
			end
			
			% update v 
			q1 = GammaG_directed(1:n, :) * (hv .* psi(:,i+1));
			q2 = GammaB_undirected(1:n, :) * (hv .* (1-sqrt_term_psi));
			v(:, i+1) = ones(n,1) - ((4*S) \ ((QL - q1 - q2) ./ v(:,i)));
			
			% compute iteration mismatch
			mismatch(i+1) = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
									GammaG_undirected, GammaB_directed, ...
									GammaB_undirected, v(:,i+1), psi(:,i+1), Vstar, A, Rt, hv);
			i = i+1;
		end
	end
	
	%% gather output data
	result = ([]);
	if mismatch(i) <= tol && ~psi_norm_big
		result.success = 1;
		result.iterations = i;
		result.mismatch = mismatch(1:i);
		result.v = v(:, i);			result.psi = psi(:, i);

		% recover VM, VA solutions (in original indexing) from v, psi
		VMsol = [v(:, i) .* Vstar_L; VG];
		result.VMsol = VMsol(idx_LG2og);
		
		Va_sol = rad2deg(A' \ asin(psi(:, i)));
		Va_sol = Va_sol - Va_sol(ref_idx);   % make sure ref bus angle = 0
		result.VAsol = Va_sol(idx_LG2og);
	else
		result.success = 0;
		if psi_norm_big 
			result.psi_norm_big = 1;
		else
			result.psi_norm_big = 0;
		end
		result.mismatch = mismatch;
	end
end

%% Local helper functions
function [mis] = compute_mismatch(P, QL, Gdiag, Bdiag_LL, GammaG_directed, ...
								  GammaG_undirected, GammaB_directed, ...
								  GammaB_undirected, v, psi, Vstar, A, Rt, hv)
	n = length(v);
	m = size(A, 1) - n;
	sqrt_term_psi = sqrt(1 - psi.^2);
	mis_P = Rt*(P - (Gdiag * (Vstar .* g(v,m)) .^2 + ...
			         GammaG_undirected * (hv .* sqrt_term_psi) + ...
				     GammaB_directed * (hv .* psi)));
    mis_Q = QL - (-Bdiag_LL * (Vstar(1:n) .* v) .^2 + ...
				   GammaG_directed(1:n, :) * (hv .* psi) - ...
				   GammaB_undirected(1:n, :) * (hv .* sqrt_term_psi));
	
	mis = norm([mis_P; mis_Q], inf);						
end

function gv = g(v, m)
	gv = [v; ones(m,1)];
end

function hv = h(v, A)
	m = size(A, 1) - length(v);
	A_from = zeros(size(A));    A_to   = zeros(size(A));
	A_from(A == 1) = 1;			A_to(A == -1) = 1;
	hv = (A_from' * g(v,m)) .* (A_to' * g(v,m));
end
