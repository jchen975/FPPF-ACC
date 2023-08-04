function [Vstar_L] = compute_open_circuit_voltage(B, VG)
	m = length(VG);
	n = size(B,1) - m;
	BLL = B(1:n, 1:n);
	BLG = B(1:n, n+1:end);
	Vstar_L = -BLL \ (BLG * VG);
end