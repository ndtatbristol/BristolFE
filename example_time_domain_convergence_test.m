%Example time-domain simulation of elastic waves using BristolFE
restoredefaultpath;
clear; 
close all;

%Material properties (SI units used throughout)
long_vel = 6000;
shear_vel = 3000;
density = 3700;

%mesh details
stress_state = 'plane strain';
% stress_state = 'plane stress';
structured_mesh = 0;
use_diagonal_lumped_mass_matrix = 0;
dof_per_node = 2; %DOF per node - usually this will be two, but in theory you could have 3 even with a 2D model if you wanted to (e.g. for out-of-plane SH waves) 

thickness = 10e-3;
width = 1e-3;

%Details of excitation - a Hanning-windowed toneburst
centre_freq = 5e6;
number_of_cycles = 5;
forcing_point = [0, 0] * thickness; %centre of bottom of L-shape
forcing_dir = 2; %vertical forcing

%Element size determined by wavelength
elements_per_wavelength = 20;

%Safety factor for time-steps
safety_factor = 3;

%How long to run model for - e.g. how many times through the thickness will
%wave propagate
max_transits = 5; 
max_time = thickness * max_transits / long_vel;

%Field output (displacments at all nodes, but not every time step
field_output_every_n_frames = inf;

%Nodes/directions where time-history will be recorded (displacements at
%every time-step)
history_point = [0, 0.5] * thickness; %halfway through plate
history_dir = 2;

%--------------------------------------------------------------------------

%Engineering constants from velocities and density
[lambda, mu] = fn_lame_from_velocities_and_density(long_vel, shear_vel, density);
[youngs_modulus, poissons_ratio] = fn_youngs_from_lame(lambda, mu);


%Work out required element size and time step (based on wavelength of bulk
%longitudinal waves at centre frequency of excitation)
wavelength = long_vel / centre_freq;
element_size = wavelength / elements_per_wavelength;

%Work out time step
time_step = element_size / long_vel / safety_factor;

%Mesh shape with triangular elements of appropriate size
corner_points = [
    -0.5 * width, 0
    -0.5 * width, thickness
     0.5 * width, thickness
     0.5 * width, 0
    ];

if structured_mesh
    [nodes, elements] = fn_rectangular_structured_mesh(corner_points([1,3],:), element_size);
else
    addpath('Mesh2d v24');
    hdata.hmax = element_size;
    options.output = false;
    [nodes, elements] = mesh2d(corner_points, [], hdata, options);
end

%Work out excitation signal
time = [0: time_step: max_time];
ct = number_of_cycles / centre_freq / 2;
forcing_functions = sin(2 * pi * centre_freq * (time - ct)) .* ...
    (1 + cos(2 * pi * centre_freq * (time - ct) / number_of_cycles)) / 2 .* ...
    (time <= number_of_cycles / centre_freq);


[f, forcing_nodes] = fn_apply_stress_along_line(nodes, corner_points(1,:), corner_points(4,:), element_size / 10, [0,1,0]);

forcing_functions = f(forcing_nodes, 2) * forcing_functions;

% [forcing_nodes, tmp] = fn_find_nodes_on_line(nodes, corner_points(1,:), corner_points(4,:), element_size / 10);
% force_nodes_x = corner_points(1,1) + (corner_points(4,1)-corner_points(1,1)) * tmp;
forcing_directions = ones(size(forcing_nodes)) * forcing_dir;
% forcing_functions = repmat(forcing_functions, [length(forcing_nodes), 1]);

[history_nodes, tmp] = fn_find_nodes_on_line(nodes, corner_points(2,:), corner_points(3,:), element_size / 10);
history_nodes_x = corner_points(2,1) + (corner_points(3,1)-corner_points(2,1)) * tmp;
history_directions = ones(size(history_nodes)) * history_dir;

%Side nodes for BCs
LH_nodes = fn_find_nodes_on_line(nodes, corner_points(1,:), corner_points(2,:), element_size / 10);
RH_nodes = fn_find_nodes_on_line(nodes, corner_points(3,:), corner_points(4,:), element_size / 10);

%Display mesh and excitation signal
figure;
display_options.node_sets_to_plot(1).nd = forcing_nodes;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = history_nodes;
display_options.node_sets_to_plot(2).col = 'g.';
display_options.node_sets_to_plot(3).nd = [LH_nodes; RH_nodes];
display_options.node_sets_to_plot(3).col = 'c>';
fn_display_result(nodes, elements, display_options);
title('Original mesh');

figure;
plot(time, forcing_functions(1,:));
title('Excitation signal');

%An n x 1 matrix defining the material of each element (here they are all
%the material number 1, so all entries in this matrix are 1)
element_materials = ones(size(elements, 1), 1);

materials(1).density = density;
switch stress_state
    case 'plane stress'
        materials(1).stiffness_matrix = fn_isotropic_plane_stress_stiffness_matrix(youngs_modulus, poissons_ratio);
    case 'plane strain'
        materials(1).stiffness_matrix = fn_isotropic_plane_strain_stiffness_matrix(youngs_modulus, poissons_ratio);
end

%FIRST BIT OF FE CALCULATION - Build global matrices
[K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_build_global_matrices(nodes, elements, element_materials, materials);

%BCs for sides - remove rows/cols from global matrices
BC_indices= [...
    fn_nodes_and_dofs_to_indices(LH_nodes, 1, global_matrix_nodes, global_matrix_dofs); ...
    fn_nodes_and_dofs_to_indices(RH_nodes, 1, global_matrix_nodes, global_matrix_dofs)];
[K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_reduce_global_matrics(K, M, Q, global_matrix_nodes, global_matrix_dofs, BC_indices);

fprintf('Size of model: %i DOF\n', size(K, 1));

%SECOND BIT OF FE CALCULATION - Time marching bit
[history_output, field_output] = fn_explicit_dynamic_solver(K, M, global_matrix_nodes, global_matrix_dofs, time, forcing_nodes, forcing_directions, forcing_functions, history_nodes, history_directions, field_output_every_n_frames, use_diagonal_lumped_mass_matrix);

if ~isempty(history_output)
    figure;
    plot(time, history_output);
    hold on;
    h = abs(fn_hilbert(sum(history_output))) / size(history_output, 1);
    h1 = round(length(h) / 2)
    [~, i1] = max(h(1:h1));
    [~, i2] = max(h(h1+1:end));
    i2 = i2 + h1;
    v = 2*thickness / (time(i2) - time(i1));
    plot(time, h, 'r');
    xlabel('Time (s)');
    title(sprintf('Measured velocity: %.1f m/s', v));
end

if ~isempty(field_output)
    stress = Q * field_output;
    figure;
    cmax = max(abs(stress(:)));
    for i = 1:size(stress,2)
        display_options.element_edge_color = 'None';
        display_options.scale_factor = 1;
        sigma_11 = stress(1:3:end, i);
        sigma_22 = stress(2:3:end, i);
        element_colour = sqrt(sigma_11 .^ 2 + sigma_22 .^ 2);
        cla;
        fn_display_result(nodes, elements, display_options, zeros(size(nodes)), element_colour)
        colorbar;
        caxis([0, cmax]);
        title('Stress field');
        pause(0.01);
    end
end


