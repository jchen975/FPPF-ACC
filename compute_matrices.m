function [GammaG_directed, GammaG_undirected, GammaB_directed, GammaB_undirected, S] ...
	= compute_matrices(Yft, Ytf, Y, Vstar_L, V_G, A)
%%
%	INPUTS: 
%		1. Yft
%		2. Ytf
%		3. Y
%		4. Vstar_L
%		5. V_G
%		6. A
%   OUTPUTS:
% 		1. GammaG_directed
%		2. GammaG_undirected
%		3. GammaB_directed
%		4. GammaB_undirected
%		5. S
%
%%
	n = length(Vstar_L);  % number of PQ buses
	
	% get the real and imaginary off-diagonal entries of Y, as arrays
	Gft = real(Yft); Bft = imag(Yft);
	Gtf = real(Ytf); Btf = imag(Ytf);
	
	% get the LL submatrix of the susceptance matrix
	BLL  = imag(Y(1:n, 1:n));
	
	% get from/to partition of the incidence matrix
	A_from = zeros(size(A)); 
	A_to   = zeros(size(A));
	A_from(A == 1) = 1;
	A_to(A == -1) = 1;
	
	% open circuit voltage magnitude for ALL buses
	Vstar = [Vstar_L; V_G];  
	
	% the Edge x Edge [Vi*Vj*G], [Vi*Vj*B] matrices for both ij and ji 
	DG_ft = sparse(diag(Gft)) * diag(A_from' * Vstar) * diag(A_to' * Vstar);
	DG_tf = sparse(diag(Gtf)) * diag(A_from' * Vstar) * diag(A_to' * Vstar);
	DB_ft = sparse(diag(Bft)) * diag(A_from' * Vstar) * diag(A_to' * Vstar);
	DB_tf = sparse(diag(Btf)) * diag(A_from' * Vstar) * diag(A_to' * Vstar);
	
	% the asymmetrically weighted incidence matrices, weighted by the diagonal
	% matrices defined above
	GammaG_directed   = sparse(A_from * DG_ft - A_to * DG_tf);
	GammaG_undirected = sparse(A_from * DG_ft + A_to * DG_tf);
	GammaB_directed   = sparse(A_from * DB_ft - A_to * DB_tf);
	GammaB_undirected = sparse(A_from * DB_ft + A_to * DB_tf);
	
	% stiffness matrix
	S = (1/4) * sparse(diag(Vstar_L)) * BLL * sparse(diag(Vstar_L));
end