function [history_output, field_output] = fn_explicit_dynamic_solver(K, M, global_matrix_nodes, global_matrix_dofs, time, forcing_nodes, forcing_dofs, forcing_functions, history_nodes, history_dofs, field_output_every_n_frames, use_diagonal_lumped_mass_matrix, varargin)
%SUMMARY
%   Solves explicit dynamic FE problem given applied displacements or
%   applied forces
%INPUTS
%   K - m x m global stiffness matrix
%   M - m x m global mass matrix
%   global_matrix_nodes, global_matrix_dofs - m element vectors decribing
%   nodes and dof represented by rows/cols in global matrices
%   time = n-element vector of times
%   forcing_nodes - p-element vector of nodes at which forcing functions
%   will be applied
%   forcing_dofs - p-element vector of directions in which forcing functions
%   will be applied
%   forcing_functions - p x n matrix of forces to apply
%   output every field_output_every_n_frames (set to inf to get no field
%   output)
%   history_nodes - q-element vector of node indices for complete 
%   time-history outputs
%   history_dofs - q-element vector of directions for complete 
%   time-history outputs
%   field_output_every_n_frames - complete displacement field will be
%   use_diagonal_lumped_mass_matrix - forces use of diagonal mass matrix,
%   which is faster, possibly less accurate
%OUTPUTS
%   field_output - m x floor(n / field_output_every_n_frames) matrix of displacements at all nodes
%   history_output - q x n matrix of time histories

%--------------------------------------------------------------------------

if isempty(varargin)
	use_gpu_if_present = 1;
else
	use_gpu_if_present = varargin{1};
end

gpu_present = fn_test_if_gpu_present_and_working;
if use_gpu_if_present && gpu_present
	use_gpu = 1;
else
	use_gpu = 0;
end

%Error checks
if size(K,1) ~= size(K, 2)
    error('K must be square matrix');
end
if size(M,1) ~= size(M, 2)
    error('M must be square matrix');
end
% if size(applied_forces, 2) ~= 2 | size(applied_displacements, 2) ~= 2
%     error('Second dimension of applied_forces and applied_displacements must be 2');
% end
% if size(K,1) ~= size(applied_forces,1) * size(applied_forces,1) | size(K,1) ~= size(applied_displacements,1) * size(applied_displacements,2)
%     error('Product of first two dimensions of applied_forces and applied_displacements matrices must equal size of K');
% end

%Convert nodes and directions to indices
forcing_indices = fn_nodes_and_dofs_to_indices(forcing_nodes, forcing_dofs, global_matrix_nodes, global_matrix_dofs)
history_indices = fn_nodes_and_dofs_to_indices(history_nodes, history_dofs, global_matrix_nodes, global_matrix_dofs)

%initialise history and field output variables
if ~history_indices
    history_output = zeros(length(history_indices), length(time));
end
if ~isinf(field_output_every_n_frames)
    field_output_ti = 1:field_output_every_n_frames:length(time);
    field_output = zeros(size(K, 1), length(field_output_ti));
else
    field_output_ti = [];
    field_output = [];
end

if use_diagonal_lumped_mass_matrix
    tmp = sum(M);
    inv_M = spdiags(1 ./ tmp.', 0, size(K,1), size(K,2));
end

u_previous = zeros(size(K, 1), 1);
u_dot_previous = zeros(size(K, 1), 1);
use_gpu
if use_gpu
	K = gpuArray(K);
	u_previous = gpuArray(u_previous);
	u_dot_previous = gpuArray(u_dot_previous);
	if use_diagonal_lumped_mass_matrix
		inv_M = gpuArray(inv_M);
	else
		M = gpuArray(M);
	end
	if ~history_indices
		history_indices = gpuArray(history_indices);
		history_output = gpuArray(history_output);
	end
	if ~isinf(field_output_every_n_frames)
		field_output = gpuArray(field_output);
	end
end

%Main time marching loop
time_step = time(2) - time(1);
t1 = clock;
for ti = 1:length(time)
    %set force at forcing node equal to excitation signal at this instant
    %in time
	if use_gpu
		f = gpuArray.zeros(size(K, 1), 1);
	else
		f = zeros(size(K, 1), 1);
	end
    if ti <= size(forcing_functions, 2)
        f(forcing_indices) = forcing_functions(:, ti);
    end
    %work out acceleration
    if use_diagonal_lumped_mass_matrix
        u_dot_dot_previous = inv_M * (f - K * u_previous);
    else
        u_dot_dot_previous = M \ (f - K * u_previous);
    end
    %work out velocity
    u_dot = u_dot_previous + time_step * u_dot_dot_previous;
    %work out displacement at next time step
    u = u_previous + time_step * u_dot;
    
    %history output
    history_output(:, ti) = u(history_indices);
    
    %field output
    [tmp, fi] = ismember(ti, field_output_ti);
    if tmp
        field_output(:, fi) = u;
    end
    
    %overwrite previous values with current ones ready for next loop
    u_previous = u;
    u_dot_previous = u_dot;
    
    %Show how far through calculation is
    fprintf('Time step %i of %i\n', ti, length(time));
end
if use_gpu
	if ~history_indices
		history_output = gather(history_output);
	end
	if ~isinf(field_output_every_n_frames)
		field_output = gather(field_output);
	end
end

t2 = etime(clock, t1);
fprintf('    ... completed in %.2f secs\n', t2);


end