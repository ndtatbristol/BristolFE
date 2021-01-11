function [f, node_list] = fn_apply_stress_along_line(nodes, p1, p2, tol, stress)
%SUMMARY
%   Converts specified constant stress to be applied along a line into
%   nodal forces
%INPUTS
%   nodes - n x 2 matrix of nodal coordinates
%   p1 - 2 element vector of coordinates of one end of line
%   p2 - 2 element vector of coordinates of other end of line
%   tol - distance specifying tolerance of distance to line for nodes
%   stress - 3 element vector stresses to apply (sigma_11, sigma_22, sigma_12)
%OUTPUTS
%   f - n x 2 matrix of nodal forces
%   node_list - list of node indices where forces were applied (useful for
%   visualisation)

%--------------------------------------------------------------------------

%error checks
if size(nodes, 2) ~= 2
    error('Nodes input must be two column matrix');
end
if length(p1(:)) ~= 2 | length(p2(:)) ~= 2
    error('Points defining line must be 2 element vectors');
end
if ~isscalar(tol)
    error('Tolerance must be scalar value');
end
if length(stress(:)) ~= 3
    error('Stress vector must have 3 elements');
end

%find nodes on line
[node_list, s] = fn_find_nodes_on_line(nodes, p1, p2, tol);

%work out vector of node spacings
[dummy, ii] = sort(s);
node_list = node_list(ii);
dxy = abs(nodes(node_list(1:end-1), :) - nodes(node_list(2:end), :));

%forces on each gap between nodes in x and y directions
fxy_gap = [dxy(:,2) * stress(1) + dxy(:,1) * stress(3), dxy(:,1) * stress(2) + dxy(:,2) * stress(3)];

%share forces on nodes at each of gaps x and y
fxy = ([fxy_gap; 0, 0] + [0, 0; fxy_gap]) / 2;

%add onto into force vector
f = zeros(size(nodes));
for ii = 1:size(fxy,1)
    f(node_list(ii), :) = fxy(ii, :);
end

end