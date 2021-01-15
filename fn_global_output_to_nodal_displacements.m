function u = fn_global_output_to_nodal_displacements(nodes, field_output, global_matrix_nodes, global_matrix_dofs)
%SUMMARY
%   Converts displacment vector (or matrix of displacement vectors) output
%   from FE calculation back to nodal displacements

u = zeros([size(nodes), size(field_output, 2)]);
for i = 1:length(global_matrix_nodes)
    j = global_matrix_nodes(i);
    k = global_matrix_dofs(i);
    u(j,k,:) = field_output(i, :);
end