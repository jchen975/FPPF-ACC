function [P, Q] = compute_power_injection(bus, gen, baseMVA, n, m)
	PD_load = bus(1:n, 3);       % real power demand at the load buses
	QD_load = bus(1:n, 4);		 % reactive ...
	PD_gen  = bus(n+1:end, 3);   % real power demand at the generator buses
	QD_gen  = bus(n+1:end, 4);   % reactive ...
	
	% if there are multiple generators connected to a single bus, then aggregate
	% the PG, QG to the first generator and set the rest to be 0 
	i = 1;   k = 1;   counter = 0;
	gen_idx = zeros(m,1); 
	while i <= m
		j = i + 1;
		if j <= m
			if gen(i, 1) == gen(j, 1) 
				multiple_gen = true;
				while multiple_gen
					if (gen(i, 1) ~= gen(j, 1)) || j > m
						multiple_gen = false;
					else
						gen(i, 2) = gen(i, 2) + gen(j, 2);     gen(j, 2) = 0.0;
						gen(i, 3) = gen(i, 3) + gen(j, 3);     gen(j, 3) = 0.0;
						counter = counter + 1;
						fprintf('Merged generator %3d PG, QG at bus %4d\n', counter, gen(i,1));
						j = j + 1;
					end
				end
			end
		end
		gen_idx(k) = i;
		k = k + 1;
		i = j;
	end
	PG_gen  = gen(gen_idx, 2);     % real power generation at the generator buses
	QG_gen  = gen(gen_idx, 3);	   % reactive...

	P  = [-PD_load; (PG_gen - PD_gen)] ./ baseMVA;
	Q  = [-QD_load; (QG_gen - QD_gen)] ./ baseMVA;
end