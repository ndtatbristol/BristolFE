function [nodes, elements] = fn_isometric_structured_mesh(corner_nodes, element_size)
%SUMMARY
%   Utility function for generating a isometric structured mesh of triangular
%   elements. Mesh will have jagged left and right edges
%INPUTS
%   corner_nodes - 2 x 2 matrix of nodal coordinates - each row is
%   coordinate of diagonally-opposite corners

%--------------------------------------------------------------------------
%first generate a rectangular one height

sin60 = sind(60);

block_size_x = abs(corner_nodes(2, 1) - corner_nodes(1, 1)) + element_size;
block_size_y = (abs(corner_nodes(2, 2) - corner_nodes(1, 2)));

%Work out how many nodes are needed in x and y
nodes_in_x_direction = ceil(block_size_x / element_size) + 1;
nodes_in_y_direction = ceil(block_size_y / element_size / sin60) + 1;

%Work out nodal coordinates
x = linspace(min(corner_nodes(:, 1)) - element_size / 2, max(corner_nodes(:, 1)) + element_size / 2, nodes_in_x_direction);
y = linspace(min(corner_nodes(:, 2)), max(corner_nodes(:, 2)), nodes_in_y_direction);
[node_x_positions, node_y_positions] = meshgrid(x, y);

%Now shuffle rows of x positions back/forward by half an element
node_x_positions(1:2:end, :) = node_x_positions(1:2:end, :) + element_size / 4;
node_x_positions(2:2:end, :) = node_x_positions(2:2:end, :) - element_size / 4;

%Work out node numbers associated with each element (a bit fiddly as you can see) 
node_numbers = reshape([1:numel(node_x_positions)], nodes_in_y_direction, nodes_in_x_direction);

element_node1a = node_numbers(1:2:end-1, 1:end-1);
element_node2a = node_numbers(2:2:end, 2:end);
element_node3a = node_numbers(2:2:end, 1:end-1);

element_node1b = node_numbers(1:2:end-1, 1:end-1);
element_node2b = node_numbers(1:2:end-1, 2:end);
element_node3b = node_numbers(2:2:end, 2:end);

element_node1c = node_numbers(2:2:end-1, 2:end);
element_node2c = node_numbers(3:2:end, 2:end);
element_node3c = node_numbers(3:2:end, 1:end-1);

element_node1d = node_numbers(2:2:end-1, 1:end-1);
element_node2d = node_numbers(2:2:end-1, 2:end);
element_node3d = node_numbers(3:2:end, 1:end-1);

%Final m x 2 matrix of x and y coordinates for each node
nodes = [node_x_positions(:), node_y_positions(:)];
%Final n x 3 matrix of 3 node numbers for each element
elements = [
    element_node1a(:), element_node2a(:), element_node3a(:)
    element_node1b(:), element_node2b(:), element_node3b(:)
    element_node1c(:), element_node2c(:), element_node3c(:)
    element_node1d(:), element_node2d(:), element_node3d(:)
    ];

end