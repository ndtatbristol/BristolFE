function indices = fn_nodes_and_dofs_to_indices(nodes, dofs, global_matrix_nodes, global_matrix_dofs)
%SUMMARY
%   Utility function for converting node/dof pairs to indices into global
%   matrices
%INPUTS
%   nodes - n x 1 element vector
%   dofs - n x 1 element vector or scalar
%   global_matrix_nodes, global_matrix_dofs - m element vectors decribing
%   nodes and dof represented by rows/cols in global matrices
%OUTPUTS
%   indices - n x 1 element vector

if isscalar(dofs)
    dofs = ones(size(nodes)) * dofs;
end

indices = zeros(length(nodes), 1);
for i = 1:length(nodes)
    tmp = find(global_matrix_nodes == nodes(i) & global_matrix_dofs == dofs(i));
    if isempty(tmp)
        warning('Node/DOF not found in global matrix');
        indices(i) = NaN;
    else
        indices(i) = tmp;
    end
end
return