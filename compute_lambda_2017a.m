function Multiplier = compute_lambda_2017a(mpc)
% This code is originally written by John W. Simpson-Porco.

%Set Continuiaton Options

define_constants;

mpc_opt = mpoption('out.all', 0, 'verbose', 0);
mpc_opt = mpoption(mpc_opt, 'cpf.stop_at', 'NOSE', 'cpf.adapt_step', 1);

%Initialize Two Cases, one base and one target
mpc_base = mpc;
mpc_target = mpc;

mpc_base.gen(:, [PG, QG]) = mpc.gen(:, [PG, QG]);
mpc_base.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]);

mpc_target.gen(:, PG) = 2*mpc.gen(:, PG);  % increased generation
mpc_target.bus(:, PD) = 2*mpc.bus(:, PD);  % and increased load
mpc_target.bus(:, QD) = 2*mpc.bus(:, QD);  % and increased load

results = runcpf(mpc_base, mpc_target, mpc_opt);

Multiplier = (1+results.cpf.max_lam);

end