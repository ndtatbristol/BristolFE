%Example time-domain simulation of elastic waves using BristolFE
%See example_time_domain_convergence_test.m for more detailed comments
restoredefaultpath;
clear; 
close all;

%Material properties (SI units used throughout)
vel = 1500;
density = 1000;

use_diagonal_lumped_mass_matrix = 1;

%Define shape of a 2D structure - a section through a plate of specified
%thickness with width 2*thickness
thickness = 5e-3;
corner_points = [
    -1, 0
    -1, 1
    1, 1
    1, 0
    ] * thickness;

%Details of excitation - a Hanning-windowed toneburst
centre_freq = 5e6;
number_of_cycles = 5;
forcing_point = [0, 0] * thickness; %centre of bottom of model

%Element size determined by wavelength
elements_per_wavelength = 10;

%Safety factor for time-steps
safety_factor = 3;

%How long to run model for - e.g. 20 x period of excitation
max_time = 1 / centre_freq * 30;

%Field output (displacements at all nodes, but not every time step - use 
%field_output_every_n_frames = inf to prevent field output).
field_output_every_n_frames = inf;

%Nodes/directions where time-history will be recorded (displacements at
%every time-step)
history_point = [0, 0.5] * thickness; %halfway through plate

%END OF INPUTS
%--------------------------------------------------------------------------

%Engineering constants from velocities and density
shear_modulus = vel ^ 2 * density;

%Work out required element size and time step (based on wavelength of bulk
%longitudinal waves at centre frequency of excitation)
% wave_vel = sqrt(youngs_modulus/density * (1-poissons_ratio) / (1+poissons_ratio) / (1-2*poissons_ratio)); %Textbook equation for bulk longitudinal wave speed
wavelength = vel / centre_freq;
element_size = wavelength / elements_per_wavelength;

%Work out time step
time_step = element_size / vel / safety_factor;

%Mesh shape with triangular elements of appropriate size
[nodes, elements] = fn_isometric_structured_mesh(corner_points([1,3],:), element_size);

%Work out excitation signal
time = [0: time_step: max_time];
ct = number_of_cycles / centre_freq / 2;
forcing_functions = sin(2 * pi * centre_freq * (time - ct)) .* ...
    (1 + cos(2 * pi * centre_freq * (time - ct) / number_of_cycles)) / 2 .* ...
    (time <= number_of_cycles / centre_freq);
forcing_nodes = fn_find_node_at_point(nodes, forcing_point, inf);

history_nodes = fn_find_node_at_point(nodes, history_point, inf);

%Display mesh
figure;
display_options.node_sets_to_plot(1).nd = forcing_nodes;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = history_nodes;
display_options.node_sets_to_plot(2).col = 'g.';
fn_display_result(nodes, elements, display_options);
title('Original mesh');

%Display excitation signal
figure;
plot(time, forcing_functions);
title('Excitation signal');

%An n x 1 matrix defining the material of each element (here they are all
%the material number 1, so all entries in this matrix are 1)
element_materials = ones(size(elements, 1), 1);

materials(1).density = density;
materials(1).stiffness_matrix = shear_modulus;

%FIRST BIT OF FE CALCULATION - Build global matrices
[K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_build_global_matrices(nodes, elements, element_materials, materials, '3');

fprintf('Size of model: %i DOF\n', size(K, 1));

%SECOND BIT OF FE CALCULATION - Time marching bit
forcing_dofs = ones(size(forcing_nodes)) * 3; %everything is in 3-direction in scalar model
history_dofs = ones(size(history_nodes)) * 3;
[history_output, field_output] = fn_explicit_dynamic_solver(K, M, global_matrix_nodes, global_matrix_dofs, time, forcing_nodes, forcing_dofs, forcing_functions, history_nodes, history_dofs, field_output_every_n_frames, use_diagonal_lumped_mass_matrix);

%Display history output
if ~isempty(history_output)
    figure;
    plot(time, history_output);
    hold on;
    h = abs(fn_hilbert(sum(history_output,1))) / size(history_output, 1);
    plot(time, h, 'r');
    xlabel('Time (s)');
end

%Animate field output
if ~isempty(field_output)
    stress = Q * field_output;
    figure;
    cmax = max(abs(stress(:)));
    for i = 1:size(stress,2)
        display_options.element_edge_color = 'None';
        display_options.scale_factor = 1;
        sigma_11 = stress(1:2:end, i);
        sigma_22 = stress(2:2:end, i);
        element_colour = sqrt(sigma_11 .^ 2 + sigma_22 .^ 2);
        cla;
        fn_display_result(nodes, elements, display_options, zeros(size(nodes)), element_colour)
        colorbar;
        caxis([0, cmax]);
        title('Stress field');
        pause(0.01);
    end
end


