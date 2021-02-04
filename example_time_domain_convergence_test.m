%Example time-domain simulation of plane elastic waves using BristolFE 
%functions in a plane strain model with symmetry boundary conditions.
%Combination of plane strain and BCs means that model should simulate plane
%bulk longitudinal waves and can be used for assessing accuracy of velocity
%as meshing parameters etc are altered.
restoredefaultpath;
clear; 
close all;

%%%%%%%%%%% MATERIAL PROPERTIES
%Material properties are defined in terms of
%ultrasonic longitudinal and shear velocities and density (these are 
%converted to Young's modulus and Poisson's ratio for FE calculation)
long_vel = 6000;
shear_vel = 3000;
density = 3700;

%%%%%%%%%%% MODEL GEOMETRY
%Here model is a simple 2D rectangle of material defined by two
%dimensions. Model extends from -width/2 to width/2 in x direction and from
%0 to thickness in y-direction
thickness = 10e-3;
width = 1e-3;

%%%%%%%%%%% MESH DETAILS
%The stress_state describes what happens in z-direction and can be either 
%'plane strain' or 'plane stress'. For ultrasonic modelling, plane strain
%is the usual assumption
stress_state = 'plane strain'; 
%The way the mesh is generated can be either 'structured' or 'unstructured'
%Here a structured mesh is made from a grid of squares each divided into
%two right angled triangles, although alternative structured meshing
%functions could be written, e.g. based on equilateral triangles. If the
%'unstructured' option is chosen, a special meshing programme is use to
%generate a mesh of triangular elements that fill the model and follow the
%geometry, which is good for models with non-rectilinear shapes.
mesh_type = 'unstructured';

%Dynamic FE simulations require mass to be considered. There are different
%ways of computing the associated mass matrix for each element in the 
%model. The 'consistent mass matrix' method forms it in a way that is 
%consistent with the way element stiffness matrix is formed; the 'lumped
%mass matrix' method just concentrates the mass of each element at its
%nodes. This results in a diagonal mass matrix. In explicit dynamic models, 
%the mass matrix has to be inverted at each time increment and it is much 
%faster to invert a diagonal matrix, even if it is a less accurate
%approximation than a consistent mass matrix. Set
%use_diagonal_lumped_mass_matrix = 1 to use a lumped mass matrix and
%use_diagonal_lumped_mass_matrix = 0 to use consistent mass matrix.
use_diagonal_lumped_mass_matrix = 1;

%The main factor governing accuracy of dynamic simulations it the element 
%size relative to the wavelength. Therefore, the element size here is
%specified in terms of elements_per_wavelength (where the reference 
%wavelength used for the calculation is that of the longitudinal wave at
%the centre frequency of excitation. Using a larger value of 
%elements_per_wavelength will result in a more accurate but slower model.
elements_per_wavelength = 20;

%%%%%%%%%%% EXCITATION AND TIME STEPPING
%In this model, the excitation is a time-varying normal stress applied
%along the lower edge of model (i.e. sigma_yy). The time variation is a 
%windowed toneburst with a specified centre frequency and number of cycles.
centre_freq = 5e6;
number_of_cycles = 5;

%An explicit dynamic model will run through a number of time steps until a
%specified end time is reached. Here the end time is defined in terms of
%the number of times a longitudinal wave will traverse the length (y
%dimension) of the model.
max_transits = 5; 

%The time-step in an explicit dynamic model must be small enough for the 
%model to be stable. The stability criteria is related to the length of
%time it takes the fastest wave to traverse the smallest element. The time
%step used is hence calculated by estimating this time and dividing by the
%safety_factor below. The safety_factor can be reduced and the model will
%run faster (because it will use fewer larger time steps to get to the same 
%end point) but it may become unstable.
safety_factor = 3;

%%%%%%%%%%% MODEL OUTPUTS

%When the model runs, the displacement at every node and every time step is
%calculated. However, only specified values are saved and output. History
%values are displacements at all time points for specified nodes; field
%outputs are displacements at all nodes at specified time points.

%Field output (displacments at all nodes, but not every time step). Use
%this for producing animations of the waves in the model. If these are not
%required set field_output_every_n_frames to inf (i.e. never give a field
%output).
field_output_every_n_frames = inf;

%In this particular model, the history output is the y-displacement of
%all nodes on the top edge of the model.

%%%%%%%%%%% END OF INPUTS
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

switch mesh_type
    case 'structured'
        [nodes, elements] = fn_rectangular_structured_mesh(corner_points([1,3],:), element_size);
    case 'unstructured'
        addpath('Mesh2d v24');
        hdata.hmax = element_size;
        options.output = false;
        [nodes, elements] = mesh2d(corner_points, [], hdata, options);
end

%Work out excitation signal
max_time = thickness * max_transits / long_vel;
time = [0: time_step: max_time];
ct = number_of_cycles / centre_freq / 2;
forcing_functions = sin(2 * pi * centre_freq * (time - ct)) .* ...
    (1 + cos(2 * pi * centre_freq * (time - ct) / number_of_cycles)) / 2 .* ...
    (time <= number_of_cycles / centre_freq);


%Locate forcing nodes and provide the time-varying force at each one (here
%it is the same at each node, but this does not have to be the case).
forcing_dir = 2; %forcing in y direction
forcing_stress = zeros(1,3);
forcing_stress(forcing_dir) = 1; %normal stress in required forcing direction
%Follow function converts a stress [sigma_xx, sigma_yy, sigma_xy] along a
%specified line into equivalent forces, f, to be applied at nodes in model
[f, forcing_nodes] = fn_apply_stress_along_line(nodes, corner_points(1,:), corner_points(4,:), element_size / 10, forcing_stress);
forcing_functions = f(forcing_nodes, 2) * forcing_functions;
forcing_directions = ones(size(forcing_nodes)) * forcing_dir;

%Locate nodes where history output is required (y-displacement on top edge
%of model)
history_dir = 2;
[history_nodes, tmp] = fn_find_nodes_on_line(nodes, corner_points(2,:), corner_points(3,:), element_size / 10);
history_nodes_x = corner_points(2,1) + (corner_points(3,1)-corner_points(2,1)) * tmp;
history_directions = ones(size(history_nodes)) * history_dir;

%This model has symmetry Boundary Conditions (BCs) on each vertical side.
%These are implemented by removing the associated rows/cols from global
%matrices before running the dynamic model.
%Locate nodes where BCs will be applied
LH_nodes = fn_find_nodes_on_line(nodes, corner_points(1,:), corner_points(2,:), element_size / 10);
RH_nodes = fn_find_nodes_on_line(nodes, corner_points(3,:), corner_points(4,:), element_size / 10);

%Display mesh
figure;
display_options.node_sets_to_plot(1).nd = forcing_nodes;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = history_nodes;
display_options.node_sets_to_plot(2).col = 'g.';
display_options.node_sets_to_plot(3).nd = [LH_nodes; RH_nodes];
display_options.node_sets_to_plot(3).col = 'c>';
fn_display_result(nodes, elements, display_options);
title('Original mesh');

%Display excitation signal
figure;
plot(time, forcing_functions(1,:));
title('Excitation signal');
xlabel('Time (s)');

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

%Apply BCs for sides - remove rows/cols from global matrices
BC_indices= [...
    fn_nodes_and_dofs_to_indices(LH_nodes, 1, global_matrix_nodes, global_matrix_dofs); ...
    fn_nodes_and_dofs_to_indices(RH_nodes, 1, global_matrix_nodes, global_matrix_dofs)];
[K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_reduce_global_matrics(K, M, Q, global_matrix_nodes, global_matrix_dofs, BC_indices);

fprintf('Size of model: %i DOF\n', size(K, 1));

%SECOND BIT OF FE CALCULATION - Time marching bit
[history_output, field_output] = fn_explicit_dynamic_solver(K, M, global_matrix_nodes, global_matrix_dofs, time, forcing_nodes, forcing_directions, forcing_functions, history_nodes, history_directions, field_output_every_n_frames, use_diagonal_lumped_mass_matrix);

%Show history output and estimate modelled wave velocity based on mean of 
%displacements at all history nodes (i.e. the mean y-displacment on the top
%surface of model).
if ~isempty(history_output)
    figure;
    plot(time, history_output);
    hold on;
    h = abs(fn_hilbert(mean(history_output)));
    plot(time, h, 'r');
    xlabel('Time (s)');
    
    %Estimate velocity by identifying difference in arrival time of waves
    %in first and second half of signal.
    h1 = round(length(h) / 2);
    [~, i1] = max(h(1:h1));
    [~, i2] = max(h(h1+1:end));
    i2 = i2 + h1;
    v = 2 * thickness / (time(i2) - time(i1));
    title(sprintf('Measured velocity: %.1f m/s', v));
end

%Animate field output
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
        fn_display_result(nodes, elements, display_options, zeros(size(nodes)), element_colour);
        colorbar;
        caxis([0, cmax]);
        title('Stress field');
        pause(0.01);
    end
end


