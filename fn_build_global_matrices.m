function [K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_build_global_matrices(nodes, elements, element_materials, materials, varargin)
%SUMMARY
%   Creates global matrices from mesh definitions
%INPUTS
%   nodes - n x 2 matrix of nodal coordinates. The row number is the node
%   number; columns 1 and 2 are the x and y coordinates of the node.
%   elements - m x 3 matrix of element nodes. The row number is the element
%   number; columns 1, 2 and 3 are the node numbers of the 3 nodes for each
%   triangular element
%   element_materials - m x 1 vector of element materials (which refer to
%   materials in 'materials' parameter
%   materials - p x 1 structured variable of materials with fields
%       materials(i).name - string giving name of material
%       materials(i).density - density of material
%       materials(i).stiffness_matrix - 3x3 or 6x6 stiffness matrix of 
%       material. 6x6 size us needed for model_dof='3';
%   [model_dof - string that is either '12' (default) or '3']
%OUTPUTS
%   K, M - global (d*n) x (d*n) stiffness and mass matrices, where d is no
%   of DOF at each node
%   Q - (e*m) x (d*n) matrix to transform displacements into stresses,
%   where e is number of stress components in each element
%--------------------------------------------------------------------------


%switch depending on DOF
if isempty(varargin)
    model_dof = '12';
else
    model_dof = varargin{1};
end

%hard-coded values - DOF per node and nodes per element
switch model_dof
    case '12'
        dofs = [1;2];
        no_stress_components = 3;
    case '3'
        dofs = 3;
        no_stress_components = 2;
end
nds_per_el = 3;

%check inputs
if size(nodes, 2) ~= 2
    error('Nodes input must be two column matrix');
end
if size(elements, 2) ~= 3
    error('Elements input must be three column matrix');
end
if (max(max(elements)) > size(nodes, 1)) | (min(min(elements)) < 1)
    error('Element node number(s) refers to undefined nodes');
end
if ~isvector(element_materials)
    error('Element materials must be a vector');
end
if length(element_materials) ~= size(elements, 1)
    error('Length of element materials vector must equal number of elements');
end
if ~isvector(materials)
    error('Materials structure must be a vector');
end
if (max(element_materials) > length(materials)) | (min(element_materials) < 1)
    error('Element material(s) refers to undefined material');
end

%no of DOFs present
dof_per_node = length(dofs);
no_nds = size(nodes,1);
no_els = size(elements,1);
tmp = repmat([1:no_nds], dof_per_node, 1);
global_matrix_nodes = tmp(:);
tmp = repmat(dofs, 1, no_nds);
global_matrix_dofs = tmp(:);
total_dof = length(global_matrix_nodes);

%build all element matrices

%there should be something else in here to check element types and only do
%els of same type together
switch model_dof
    case '12'
        [K_flat, M_flat, Q_flat] = fn_get_element_matrices(nodes, elements, element_materials, materials);
    case '3'
        [K_flat, M_flat, Q_flat] = fn_get_element_matrices_scalar(nodes, elements, element_materials, materials);
end
disp('Global matrices ...');
t1 = clock;

%build the matrices
[i1, i2] = fn_global_indices_for_flat_K_and_M(elements, dof_per_node, nds_per_el);
M = sparse(i1(:), i2(:), M_flat(:), total_dof, total_dof);
K = sparse(i1(:), i2(:), K_flat(:), total_dof, total_dof);
[i1, i2] = fn_global_indices_for_flat_Q(elements, dof_per_node, nds_per_el, no_stress_components);
Q = sparse(i1(:), i2(:), Q_flat(:), no_els * no_stress_components, total_dof);
Q_matrix_el = ones(no_stress_components, 1) * [1:size(elements,1)];
Q_matrix_el = Q_matrix_el(:);
Q_matrix_stress_dir = repmat([1:no_stress_components]',size(elements,1), 1);
disp(sprintf('    ... built in %.2f secs', etime(clock, t1)));

end

function [i1, i2] = fn_global_indices_for_flat_K_and_M(elements, dof_per_node, nds_per_el)
nd1 = ones(dof_per_node,1) * [1:nds_per_el];
nd1 = nd1(:);
[nd2, nd1] = meshgrid(nd1, nd1);
nd1 = nd1(:);
nd2 = nd2(:);

dof1 = [1:dof_per_node]' * ones(1, nds_per_el);
dof1 = dof1(:);
[dof2, dof1] = meshgrid(dof1, dof1);
dof1 = dof1(:);
dof2 = dof2(:);

i1 = elements(:,nd1)';
i2 = elements(:,nd2)';

for ii = 1:size(i1,1)
    i1(ii,:) = i1(ii,:) * dof_per_node - dof_per_node + 1 + dof1(ii) - 1;
    i2(ii,:) = i2(ii,:) * dof_per_node - dof_per_node + 1 + dof2(ii) - 1;
end
end

function [i1, i2] = fn_global_indices_for_flat_Q(elements, dof_per_node, nds_per_el, no_stress_components)
nd2 = ones(dof_per_node,1) * [1:nds_per_el];
nd2 = nd2(:);

[nd2, stress_comps] = meshgrid(nd2, [1:no_stress_components]);
stress_comps = stress_comps(:);
nd2 = nd2(:);

dof2 = [1:dof_per_node]' * ones(1, nds_per_el);
dof2 = dof2(:);
[dof2, dummy] = meshgrid(dof2, [1:no_stress_components]);
dof2 = dof2(:);

i1 = ones(no_stress_components * nds_per_el * dof_per_node,1) * [1:size(elements,1)];
i2 = elements(:,nd2)';
for ii = 1:size(i1,1)
    i1(ii,:) = (i1(ii,:) - 1) * no_stress_components + stress_comps(ii);
    i2(ii,:) = (i2(ii,:) - 1) * dof_per_node + dof2(ii);
end
end

function [K, M, Q] = fn_get_element_matrices(nodes, elements, element_materials, materials)
%Cutdown version of SAFE code for UG FE - constant strain triangles, 2 DOF
%per node

disp('Element matrices ...');
t1 = clock;

rho = [materials(element_materials).density];

x = zeros(3, size(elements,1));
y = zeros(3, size(elements,1));
D = zeros(9, size(elements,1));

for ii = 1:3
    x(ii, :) = nodes(elements(:, ii), 1);
    y(ii, :) = nodes(elements(:, ii), 2);
end
for ii = 1:length(materials)
    if size(materials(ii).stiffness_matrix, 1) ~= 3 && size(materials(ii).stiffness_matrix, 2) ~= 3
        error('Material stiffness matrices must be 3x3 for elastic elements elements')
    end
    materials(ii).stiffness_matrix = materials(ii).stiffness_matrix(:);
end

D = [materials(element_materials).stiffness_matrix];

J = x(1, :) .* y(2, :)  -  x(2, :) .* y(1, :)  -  x(1, :) .* y(3, :)  +  x(3, :) .* y(1, :)  +  x(2, :) .* y(3, :)  -  x(3, :) .* y(2, :);
K = zeros(36, length(J));
K(1, :) = (((D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(3, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(7, :) = (((D(6, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(4, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(13, :) = (((D(3, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(19, :) = (((D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(6, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(4, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(25, :) = (((D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(3, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(31, :) = (((D(6, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(4, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(2, :) = (((D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(2, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(8, :) = (((D(5, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(6, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(14, :) = (((D(2, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(20, :) = (((D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(5, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(6, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(26, :) = (((D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(2, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(32, :) = (((D(5, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(6, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(3, :) = (((D(3, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2  -  (((D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2;
K(9, :) = (((D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2  -  (((D(6, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(4, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2;
K(15, :) = (((D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2  -  (((D(3, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2;
K(21, :) = (((D(6, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(4, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2  -  (((D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2;
K(27, :) = (((D(3, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2  -  (((D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2;
K(33, :) = (((D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2  -  (((D(6, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(4, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2;
K(4, :) = (((D(2, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2  -  (((D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2;
K(10, :) = (((D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2  -  (((D(5, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(6, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2;
K(16, :) = (((D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2  -  (((D(2, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2;
K(22, :) = (((D(5, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(6, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2  -  (((D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2;
K(28, :) = (((D(2, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2  -  (((D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2;
K(34, :) = (((D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2  -  (((D(5, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(6, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2;
K(5, :) = (((D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(3, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(11, :) = (((D(6, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(4, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(17, :) = (((D(3, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(23, :) = (((D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(6, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(4, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(29, :) = (((D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(3, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(35, :) = (((D(6, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(4, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(6, :) = (((D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(2, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(12, :) = (((D(5, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(6, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
K(18, :) = (((D(2, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(24, :) = (((D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(5, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(6, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
K(30, :) = (((D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(2, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
K(36, :) = (((D(5, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(6, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
M = zeros(36, length(J));
M(1, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(13, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(25, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(8, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(20, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(32, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(3, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(15, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(27, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(10, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(22, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(34, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(5, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(17, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(29, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(12, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(24, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(36, :) = (J(1, :) .* rho(1, :)) ./ 12;
Q = zeros(18, length(J));
Q(1, :) = (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)  -  (D(7, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(4, :) = (D(7, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)  -  (D(4, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(7, :) = (D(7, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(10, :) = (D(4, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(7, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(13, :) = (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)  -  (D(7, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :);
Q(16, :) = (D(7, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)  -  (D(4, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :);
Q(2, :) = (D(2, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)  -  (D(8, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(5, :) = (D(8, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)  -  (D(5, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(8, :) = (D(8, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(2, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(11, :) = (D(5, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(8, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(14, :) = (D(2, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)  -  (D(8, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :);
Q(17, :) = (D(8, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)  -  (D(5, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :);
Q(3, :) = (D(3, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)  -  (D(9, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(6, :) = (D(9, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)  -  (D(6, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(9, :) = (D(9, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(3, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(12, :) = (D(6, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(9, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(15, :) = (D(3, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)  -  (D(9, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :);
Q(18, :) = (D(9, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)  -  (D(6, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :);

disp(sprintf('    ... built in %.2f secs', etime(clock, t1)));
end

function [K, M, Q] = fn_get_element_matrices_scalar(nodes, elements, element_materials, materials)
%Cutdown version of SAFE code for UG FE - constant strain triangles, 2 DOF
%per node

disp('Element matrices ...');
t1 = clock;

rho = [materials(element_materials).density];

x = zeros(3, size(elements,1));
y = zeros(3, size(elements,1));
D = zeros(36, size(elements,1));

for ii = 1:3
    x(ii, :) = nodes(elements(:, ii), 1);
    y(ii, :) = nodes(elements(:, ii), 2);
end

for ii = 1:length(materials)
    if ~isscalar(materials(ii).stiffness_matrix)
        error('Material stiffness matrices must be 1x1 for scalar elements');
    end
%     materials(ii).stiffness_matrix = materials(ii).stiffness_matrix(:);
end

G = [materials(element_materials).stiffness_matrix];
% J = abs(x(1, :) .* y(2, :)  -  x(2, :) .* y(1, :)  -  x(1, :) .* y(3, :)  +  x(3, :) .* y(1, :)  +  x(2, :) .* y(3, :)  -  x(3, :) .* y(2, :));
% K = zeros(9, length(J));
% K(1, :) = (((D(1, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(1, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
% K(4, :) = (((D(1, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(1, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
% K(7, :) = (((D(1, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(1, :) .* (x(2, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(2, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
% K(2, :) = (((D(1, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2  -  (((D(1, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2;
% K(5, :) = (((D(1, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2  -  (((D(1, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2;
% K(8, :) = (((D(1, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2  -  (((D(1, :) .* (x(1, :)  -  x(3, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(3, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2;
% K(3, :) = (((D(1, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(2, :)  -  x(3, :))) ./ 2  -  (((D(1, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(2, :)  -  y(3, :))) ./ 2;
% K(6, :) = (((D(1, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(3, :))) ./ 2  -  (((D(1, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(3, :))) ./ 2;
% K(9, :) = (((D(1, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (x(1, :)  -  x(2, :))) ./ 2  -  (((D(1, :) .* (x(1, :)  -  x(2, :))) ./ J(1, :)  -  (D(1, :) .* (y(1, :)  -  y(2, :))) ./ J(1, :)) .* (y(1, :)  -  y(2, :))) ./ 2;
% M = zeros(9, length(J));
% M(1, :) = (J(1, :) .* rho(1, :)) ./ 12;
% M(4, :) = (J(1, :) .* rho(1, :)) ./ 24;
% M(7, :) = (J(1, :) .* rho(1, :)) ./ 24;
% M(2, :) = (J(1, :) .* rho(1, :)) ./ 24;
% M(5, :) = (J(1, :) .* rho(1, :)) ./ 12;
% M(8, :) = (J(1, :) .* rho(1, :)) ./ 24;
% M(3, :) = (J(1, :) .* rho(1, :)) ./ 24;
% M(6, :) = (J(1, :) .* rho(1, :)) ./ 24;
% M(9, :) = (J(1, :) .* rho(1, :)) ./ 12;
J = x(1, :) .* y(2, :)  -  x(2, :) .* y(1, :)  -  x(1, :) .* y(3, :)  +  x(3, :) .* y(1, :)  +  x(2, :) .* y(3, :)  -  x(3, :) .* y(2, :);
K = zeros(9, length(J));
K(1, :) = (G .* (x(2, :)  -  x(3, :)) .^ 2) ./ (2 .* J(1, :))  +  (G .* (y(2, :)  -  y(3, :)) .^ 2) ./ (2 .* J(1, :));
K(2, :) =  -  (G .* (x(1, :)  -  x(3, :)) .* (x(2, :)  -  x(3, :))) ./ (2 .* J(1, :))  -  (G .* (y(1, :)  -  y(3, :)) .* (y(2, :)  -  y(3, :))) ./ (2 .* J(1, :));
K(3, :) = (G .* (x(1, :)  -  x(2, :)) .* (x(2, :)  -  x(3, :))) ./ (2 .* J(1, :))  +  (G .* (y(1, :)  -  y(2, :)) .* (y(2, :)  -  y(3, :))) ./ (2 .* J(1, :));
K(4, :) =  -  (G .* (x(1, :)  -  x(3, :)) .* (x(2, :)  -  x(3, :))) ./ (2 .* J(1, :))  -  (G .* (y(1, :)  -  y(3, :)) .* (y(2, :)  -  y(3, :))) ./ (2 .* J(1, :));
K(5, :) = (G .* (x(1, :)  -  x(3, :)) .^ 2) ./ (2 .* J(1, :))  +  (G .* (y(1, :)  -  y(3, :)) .^ 2) ./ (2 .* J(1, :));
K(6, :) =  -  (G .* (x(1, :)  -  x(2, :)) .* (x(1, :)  -  x(3, :))) ./ (2 .* J(1, :))  -  (G .* (y(1, :)  -  y(2, :)) .* (y(1, :)  -  y(3, :))) ./ (2 .* J(1, :));
K(7, :) = (G .* (x(1, :)  -  x(2, :)) .* (x(2, :)  -  x(3, :))) ./ (2 .* J(1, :))  +  (G .* (y(1, :)  -  y(2, :)) .* (y(2, :)  -  y(3, :))) ./ (2 .* J(1, :));
K(8, :) =  -  (G .* (x(1, :)  -  x(2, :)) .* (x(1, :)  -  x(3, :))) ./ (2 .* J(1, :))  -  (G .* (y(1, :)  -  y(2, :)) .* (y(1, :)  -  y(3, :))) ./ (2 .* J(1, :));
K(9, :) = (G .* (x(1, :)  -  x(2, :)) .^ 2) ./ (2 .* J(1, :))  +  (G .* (y(1, :)  -  y(2, :)) .^ 2) ./ (2 .* J(1, :));
M = zeros(9, length(J));
M(1, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(2, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(3, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(4, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(5, :) = (J(1, :) .* rho(1, :)) ./ 12;
M(6, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(7, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(8, :) = (J(1, :) .* rho(1, :)) ./ 24;
M(9, :) = (J(1, :) .* rho(1, :)) ./ 12;

%these need checking. In theory [sigma_23, sigma_13] = [u1, v1, u2, v2, u3, v3] * Q
Q = zeros(6, length(J));
Q(1, :) =  - (G .* (x(2, :)  -  x(3, :))) ./ J(1, :);
Q(2, :) = (G .* (y(2, :)  -  y(3, :))) ./ J(1, :);
Q(3, :) = (G .* (x(1, :)  -  x(3, :))) ./ J(1, :);
Q(4, :) =  - (G .* (y(1, :)  -  y(3, :))) ./ J(1, :);
Q(5, :) =  - (G .* (x(1, :)  -  x(2, :))) ./ J(1, :);
Q(6, :) = (G .* (y(1, :)  -  y(2, :))) ./ J(1, :);
disp(sprintf('    ... built in %.2f secs', etime(clock, t1)));
end