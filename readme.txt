BristolFE - Paul Wilcox
=======================

This repository contains a number of functions and example scripts for performing basic static or dynamic 2D Finite Element (FE) simulations. The only element type supported is a constant strain triangle.

Files named fn_*.m are function files and in general should not be altered. Either open the file in the Matlab editor or type help [filename] to see the help associated with a function.

Files named example_*.* are example Matlab scripts that illustrate how the functions should be called to perform typical FE simulations.

The subfolder Mesh2D v24 contains 2D meshing routines written by Darren Engwirda. These are included in the repo to ensure compatibility but are available directly from Github (https://github.com/dengwirda/mesh2d) and Matlab File Exchange (https://uk.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation).

OVERVIEW
========

The core FE functions are 

[K, M, Q, global_matrix_nodes, global_matrix_dofs] = fn_build_global_matrices(nodes, elements, element_materials, materials)
This builds the global stiffness (K) and mass (M) matrices for the FE model. The matrix Q is used to convert the nodel displacments form a model into stress components [sigma_xx, sigma_yy, sigma_xy] within each element. The outputs global_matrix_nodes and global_matrix_dofs list the nodes and DOF associated with each row/col in the global matrices.

[u, f] = fn_static_solver(K, applied_forces, applied_displacements)
This performs a static analysis of a model subject to applied forces and/or displacements. The outputs are the displacement (u) and external force (f) at each node.

[history_output, field_output] = fn_explicit_dynamic_solver(K, M, global_matrix_nodes, global_matrix_dofs, time, forcing_nodes, forcing_dofs, forcing_functions, history_nodes, history_dofs, field_output_every_n_frames, use_diagonal_lumped_mass_matrix)
This performs an explicit time-marching simulation when the specified time-varying forcing_functions are applied to the specified node DOFs. The output can be time histories of displacements at specified nodes or the complete displacment field at specified time steps.

The other functions are mostly utilities to perform common tasks to support the core FE functions:

fn_display_result - use for plotting mesh and results

fn_apply_stress_along_line - works out nodal forces equivalent to specified stress
fn_find_node_at_point - finds nearest node to a point
fn_find_nodes_on_line - find nearest nodes to a line
fn_global_output_to_nodal_displacements - converts raw displacement output vector from simulation back into two-column matrix of the displacement for each node in each DOF

fn_isotropic_plane_strain_stiffness_matrix - plane strain stiffness matrix for isotropic material based on Young's modulus and Poisson's ratio
fn_isotropic_plane_stress_stiffness_matrix - plane stress stiffness matrix for isotropic material based on Young's modulus and Poisson's ratio
fn_lame_from_velocities_and_density - converts ultrasonic longitudinal and shear velocities and density into Lame elastic constants 
fn_youngs_from_lame - converts Lame elastic constants into Young's modulus and Poisson's ratio

fn_nodes_and_dofs_to_indices - converts lists of nodes and DOFs into indices into global matrices
fn_rectangular_structured_mesh - generates rectangular structured mesh made of identical right-angle-triangle-shaped elements
fn_reduce_global_matrics - removes rows/cols from global matrices for specified nodes/DOFs. Typically used for imposing displacement BCs

fn_hilbert - performs Hilbert transform of time-domain signals



