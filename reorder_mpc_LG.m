function [load_index, gen_index, gen_bus_idx, edge_index, direction_index, ...
		  bus_idx_reverse, mpc_reordered, A, Y, Yft, Ytf] = reorder_mpc_LG(mpc)
% This script performs LG-reordering of buses, by which I mean, it reorders the 
% buses of the system placing loads (PQ) first and generators (PV and slack) 
% second. It also generates a vector of indices for edges which indicates the 
% type of edge (1=PQ-PQ, 2=PV-PQ, 3=PV-PV).

% Modified in 2021 from the code written by John W. Simpson-Porco. 

%   Inputs:
%   -----------------------------------
%   1. mpc                  the mpc file that does not satisfy the LG paritioning 


%   Outputs:
%   -----------------------------------
%   1. load_index           vector of bus indices for PQ buses
%   2. gen_index            vector of bus indices for PV buses
%   3. edge_index           edge-vector with elements {1,2,3}, indicating PQ-PQ, PV-PQ, or PV-PV edge
%   5. bus_idx_reverse      permutation to recover original (internal, sequential) mpc indexing
%   4. direction_index      edge_vector with elements {1,-1}, indicating whether branch direction in MATPOWER goes with or against the labeling convention in the paper.
%   6. mpc_reordered        the mpc file with reordered buses (WARNING: bus ordering will be non-sequential)
%   6. A                    the incidence matrix, with reordered buses (edges are NOT reordered)
%   7. Y                    the admittance matrix, with reordered buses
%   8. Yft
%   9. Ytf

GEN_STATUS = 8;
mpc = ext2int(mpc);

mpc.gen(:, GEN_STATUS) = 1;						% Force all generators to be on
[Y, Yft, Ytf, A] = makeYbus_modified(mpc);		% Build Y_bus and node-edge admittance matrix
E = size(A,2);									% Number of branches in network

%% Get Initial Indices for Node Types
old_index_set = mpc.bus(:,1);                         % Get original bus indices
% gen_index = mpc.gen(:,1);                             % Get PV/slack indices
[gen_index, gen_bus_idx] = unique(mpc.gen(:,1), 'stable');
load_index = setdiff(old_index_set,gen_index);        % Get PQ indices
bus_index_reo = [load_index;gen_index];               % Construct new bus indices
[~, bus_idx_reverse] = sort(bus_index_reo);

% [REF, PV, PQ] = bustypes(mpc.bus, mpc.gen);           % Alternative code
% gen_index2 = [PV;REF];
% load_index2 = PQ;
% bus_index_reo2 = [load_index2;gen_index2];

%% Re-order Incidence, Admittance, MPC file
A = A(bus_index_reo,:);
Y = Y(bus_index_reo,bus_index_reo);
mpc_reordered = mpc;
mpc_reordered.bus = mpc.bus(bus_index_reo,:);

%% Edge Directions

% First, define two subfunctions. The funtion isPQ (resp. isPV) returns a 1 if the index 'ind' is a PQ (resp. PV or slack) bus.

isPQ = @(ind) true-isempty((find(load_index==ind)));    
isPV = @(ind) true-isempty((find(gen_index==ind)));     

%Initialize vectors

edge_index = zeros(E,1);  % 1 if PQ-PQ branch, 2 if PV-PQ branch, and 3 if PV-PV branch
direction_index = zeros(E,1);

for i=1:E
    from_to_bus_idx=mpc.branch(i,1:2);                    % Get sending/receiving end indices of current branch
    if isPQ(from_to_bus_idx(1)) && isPQ(from_to_bus_idx(2))    % If its a PQ-PQ branch
		edge_index(i) = 1;
        direction_index(i) = 1;
	elseif (isPQ(from_to_bus_idx(1)) && isPV(from_to_bus_idx(2))) || (isPV(from_to_bus_idx(1)) && isPQ(from_to_bus_idx(2))) %If its a PV-PQ branch
        edge_index(i) = 2;
        if isPQ(from_to_bus_idx(1)) && mpc.branch(i,1) == 1
			direction_index(i) = -1;
		elseif isPQ(from_to_bus_idx(2)) && mpc.branch(i,2) == 1 
			direction_index(i) = -1;
		else
			direction_index(i) = 1;
        end
    else   % If its a PV-PV branch
        edge_index(i) = 3;
        direction_index(i) = 1;
    end
end


end

