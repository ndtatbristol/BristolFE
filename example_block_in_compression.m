%Static FE validation script

%Generates a rectangular block supported on a frictionless flat surface and
%subjected to either a constant normal stress or constant displacement on
%the top surface. 

%The mesh (i.e. node positions and which nodes form corners of each
%element) in this model is generated completely manually within this script
%because shape is so simple. For anything more complicated, use proper
%meshing program (see "example_arched_bridge.m")

%Always a good idea to clear everything at the start of a script!
clear; %clear all variables
close all; %close all windows
clc; %clear the command screen

%--------------------------------------------------------------------------
%START OF USER INPUTS
%--------------------------------------------------------------------------
%Geometry and material properties are input in this section

%Material properties (SI units used throughout)
youngs_modulus = 70e9;
poissons_ratio = 1/3;
density = 2700;

%Stress state - can be either 'plane stress' or 'plane strain'
stress_state = 'plane stress';

%Block dimensions
block_length = 0.1;
block_depth = 0.1;

%Element size
element_size = 0.005;

%Applied loading / displacement (uncomment lines accordingly)
load_case = 'stress';
applied_stress = -7e8; %should give displacement of approx -block_depth / 100

% load_case = 'displacement';
% applied_displacement = -block_depth / 100; %  this is 1% strain, so should give sigma_yy = youngs_modulus * 0.01 =7e8 everywhere and zero for other stress components

%--------------------------------------------------------------------------
%END OF USER INPUTS
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%PRE-PROCESSING
%--------------------------------------------------------------------------

%Convert material properties into correct form for BristolFE 
    %Note that BristolFE allows more than one different material to be used in
    %the same model although only one is used here.
    materials(1).density = density;
    switch stress_state
        case 'plane stress'
            materials(1).stiffness_matrix = fn_isotropic_plane_stress_stiffness_matrix(youngs_modulus, poissons_ratio);
        case 'plane strain'
            materials(1).stiffness_matrix = fn_isotropic_plane_strain_stiffness_matrix(youngs_modulus, poissons_ratio);
    end

%--------------------------------------------------------------------------
%Generate the mesh to obtain (at the end of this section)
corner_nodes = [
    0, 0
    block_length, block_depth];

[nodes, elements] = fn_rectangular_structured_mesh(corner_nodes, element_size);

%--------------------------------------------------------------------------
%A few other things to be defined to complete mesh definition 

%An n x 1 matrix defining the material of each element (here they are all
%the material number 1, so all entries in this matrix are 1)
element_materials = ones(size(elements, 1), 1);

%Find nodes on frictionless surface (i.e. bottom of block)
fixed_nodes = fn_find_nodes_on_line(nodes, [0, 0], [block_length, 0], element_size / 10);

%Display mesh
figure;
display_options.node_sets_to_plot(1).nd = fixed_nodes;
display_options.node_sets_to_plot(1).col = 'g.';
fn_display_result(nodes, elements, display_options);
title('Original mesh');

%--------------------------------------------------------------------------
%MAIN FE CALCULATIONS
%--------------------------------------------------------------------------

%FIRST BIT OF FE CALCULATION - Build global matrices
[K, M, Q] = fn_build_global_matrices(nodes, elements, element_materials, materials);
fprintf('Size of model: %i DOF\n', size(K,1));

%APPLY BOUNDARY CONDITIONS AND LOADING
%Initially define all nodal displacements to be NaN (NaN indicates 
%unknown displacement) and all applied forces to be zero
applied_displacements = NaN(size(nodes));
applied_forces = zeros(size(nodes));

%Apply frictionless boundary condition to base (no displacement in vertical
%direction, i.e. for DOF = 2)
applied_displacements(fixed_nodes, 2) = 0;

%Apply loading / displacement
switch load_case 
    case 'stress'
        [applied_forces, forcing_nodes] = fn_apply_stress_along_line(nodes, [0, block_depth], [block_length, block_depth], element_size / 10, [0, applied_stress, 0]);
    case 'displacement'
        forcing_nodes = fn_find_nodes_on_line(nodes, [0, block_depth], [block_length, block_depth], element_size / 10);
        applied_displacements(forcing_nodes, 2) = applied_displacement;
end

%SECOND BIT OF FE CALCULATION - Solve for unknown displacements at free 
%nodes and unknown external forces at fixed nodes
[u, f] = fn_static_solver(K, applied_forces, applied_displacements);

%THIRD BIT OF FE CALCULATION - Get stresses in each element from
%displacements
[sigma_xx, sigma_yy, sigma_xy] = fn_stress_from_disp(Q, u);

%--------------------------------------------------------------------------
%POST-PROCESSING
%--------------------------------------------------------------------------
%Display results
max_stress = max(max([sigma_xx, sigma_yy, sigma_xy]));
min_stress = min(min([sigma_xx, sigma_yy, sigma_xy]));

figure;
display_options.element_edge_color = 'None';
display_options.node_sets_to_plot(2).nd = forcing_nodes;
display_options.node_sets_to_plot(2).col = 'r.';

subplot(1,3,1);
element_colour = sigma_xx;
fn_display_result(nodes, elements, display_options, u, element_colour);
colorbar;
caxis([min_stress, max_stress]);
title('{\sigma}_{xx}');

subplot(1,3,2);
element_colour = sigma_yy;
fn_display_result(nodes, elements, display_options, u, element_colour);
colorbar;
caxis([min_stress, max_stress]);
title('{\sigma}_{yy}');

subplot(1,3,3);
element_colour = sigma_xy;
fn_display_result(nodes, elements, display_options, u, element_colour);
colorbar;
caxis([min_stress, max_stress]);
title('{\sigma}_{xy}');


