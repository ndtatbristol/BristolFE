%Static FE script using separate mesh generator functions

%Generates an arched bridge shaped structure, built in at ends and applies
%either uniform or point load on bridge deck. Output is display of
%distorted shape (with exagerated displacements) coloured to show stress
%level and the displacement at the central point on the bridge deck.

%Always a good idea to clear everything at the start of a script!
clear; %clear all variables
close all; %close all windows
clc; %clear the command screen

%--------------------------------------------------------------------------
%START OF USER INPUTS
%--------------------------------------------------------------------------
%Geometry and material properties are input in this section

%material properties (SI units used throughout)
youngs_modulus = 210e9;
poissons_ratio = 0.3;
density = 8900;

%Stress state - can be either 'plane stress' or 'plane strain'
stress_state = 'plane stress';

%Bridge dimensions
bridge_span = 10;
bridge_max_depth = 4;
bridge_min_depth = 1;

%Element size
element_size = 0.2;

%Loading
load_case = 'distributed'; %can be 'distributed' or 'point' in this example
total_load_per_unit_thickness = -1e6; %negative because it is downwards

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
%Generate mesh from geometry data

%At end of this section there need to be two matrices:
%   "nodes": an m x 2 matrix of (x,y) coordinates for the m nodes
%   "elements": an n x 3 matrix of nodes at 3 corners of n triangular elements

%First need to work out coordinates of the border of a polygon that defines
%bridge shape - this will be input to meshing program
    %A bit of trig to needed work out where centre of (circular) arch is and
    %what angular range it is over based on supplied geometric data
    d = bridge_max_depth - bridge_min_depth;
    r = (d ^ 2 + bridge_span ^ 2 / 4) / (2 * d);
    yc = -bridge_min_depth - r;
    theta_max = asin(bridge_span / 2 / r);
    arc_len = 2 * theta_max * r;
    theta_pts = ceil(arc_len / element_size);
    theta = linspace(-theta_max, theta_max, theta_pts)';
    %Calculate coordiates of points defining arch
    x = r * sin(theta);
    y = yc + r * cos(theta);
    %Add two more points to define ends of top of bridge
    x = [x; bridge_span / 2; -bridge_span / 2];
    y = [y; 0; 0];
%We now have coordinates (x and y) of polygon defining bridge's shape.

%Add directory of meshing functions to Matlab path and call the meshing
%program to do the work and create a mesh inside the polygon defining the
%bridge's shape

addpath('Mesh2d v24');
hdata.hmax = element_size;
options.output = false;
[nodes, elements] = mesh2d([x, y], [], hdata, options);

%--------------------------------------------------------------------------
%A few other things to be defined to complete model definition 

%An n x 1 matrix defining the material of each element (here they are all
%the material number 1, so all entries in this matrix are 1)
element_materials = ones(size(elements, 1), 1);

%Find nodes at either end of bridge, which will be fixed to represent
%built-in boundary conditions
fixed_nodes_left = fn_find_nodes_on_line(nodes, [-bridge_span / 2, 0], [-bridge_span / 2, -bridge_max_depth], element_size / 10);
fixed_nodes_right = fn_find_nodes_on_line(nodes, [bridge_span / 2, 0], [bridge_span / 2, -bridge_max_depth], element_size / 10);
fixed_nodes = [fixed_nodes_left; fixed_nodes_right];

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
fprintf('Size of model: %i DOF\n', size(K, 1));

%APPLY BOUNDARY CONDITIONS AND LOADING
%Initially define all nodal displacements to be NaN (NaN indicates 
%unknown displacement) and all applied forces to be zero
applied_displacements = NaN(size(nodes));
applied_forces = zeros(size(nodes));

%Apply built-in boundary condition (zero displacements) by specifying that 
%applied-displacements are zero at fixed nodes
applied_displacements(fixed_nodes, :) = 0;

%Find number of node at mid-point of span (needed for later for output
%anyway, even if load not applied at this point)
centre_node = fn_find_node_at_point(nodes, [0, 0], inf);

%Apply appropriate load to top deck
switch load_case
    case 'point'
        applied_forces(centre_node, 2) = total_load_per_unit_thickness;
        forcing_nodes = centre_node;
    case 'distributed'
        stress_per_unit_thickness = total_load_per_unit_thickness / bridge_span;
        [applied_forces, forcing_nodes] = fn_apply_stress_along_line(nodes, [-bridge_span / 2, 0], [bridge_span / 2, 0], element_size / 10, [0, stress_per_unit_thickness, 0]);
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

figure;
display_options.element_edge_color = 'None';
display_options.node_sets_to_plot(2).nd = forcing_nodes;
display_options.node_sets_to_plot(2).col = 'r.';
display_options.scale_factor = 400;
element_colour = sqrt(sigma_xx .^ 2 + sigma_yy .^ 2);
fn_display_result(nodes, elements, display_options, u, element_colour)
colorbar;
title('{({\sigma}_{xx}^2+{\sigma}_{yy}^2)}^{0.5}');

%Work out vertical displacement at mid-point of deck
centre_displacement = u(centre_node, 2);
fprintf('Centre displacement: %.2f mm\n', centre_displacement * 1e3); %Note conversion to mm is only performed for displat purposes