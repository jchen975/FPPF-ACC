function [A_L_ll_plus,A_L_ll_minus,A_L_gl_minus,egg,edge_index_reverse] ... 
		= incidence_matrix_partitions(n, edge_index, A)

% This file operates on the LG-reordered incidence matrix.

ell_index = find(edge_index==1);
egl_index = find(edge_index==2);
egg_index = find(edge_index==3);

% ell = length(ell_index);
% egl = length(egl_index);
egg = length(egg_index);
% E = ell + egl + egg;

A_L = A(1:n,:);
% A_G = A(n+1:n+m,:);

A_L_ll = A_L(:,ell_index);
A_L_gl = A_L(:,egl_index);
% A_L_gg = A_L(:,egg_index);
% A_G_ll = A_G(:,ell_index);
% A_G_gl = A_G(:,egl_index);
% A_G_gg = A_G(:,egg_index);

edge_index_reo = [ell_index;egl_index;egg_index];
[~,edge_index_reverse] = sort(edge_index_reo);
%A_total_reo = [A_L_ll,A_L_gl,A_L_gg;A_G_ll,A_G_gl,A_G_gg];

% AbsA_L = abs(A_L);
A_L_ll_plus = double(A_L_ll>0);
A_L_ll_minus = double(A_L_ll<0);
A_L_gl_minus = abs(A_L_gl);

end