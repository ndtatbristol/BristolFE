function [K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_reduce_global_matrics(K, M, Q, global_matrix_nodes, global_matrix_dofs, indices_to_remove)
%SUMMARY
%   Remove rows/cols from global matrices

K(indices_to_remove, :) = [];
K(:, indices_to_remove) = [];
M(indices_to_remove, :) = [];
M(:, indices_to_remove) = [];
Q(:, indices_to_remove) = [];
global_matrix_nodes(indices_to_remove) = [];
global_matrix_dofs(indices_to_remove) = [];

end