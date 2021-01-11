function [u, f] = fn_static_solver(K, applied_forces, applied_displacements)
%SUMMARY
%   Solves linear static FE problem given applied displacements and
%   applied forces
%INPUTS
%   K - m x m global stiffness matrix
%   applied_forces - m/2 x 2 matrix of forces to apply at nodes
%   and DOF in model; zeros mean no external force is applied
%   applied_displacements - m/2 x 2 matrix of displacements to apply at
%   corresponding nodes and DOF in model; zeros mean displacement is
%   constrained to zero, use NaN for unconstrained displacements
%OUTPUTS
%   u - m/2 x 2 matrix of displacements at all nodes
%   f - m/2 x 2 matrix of external forces at all nodes (e.g. reaction forces
%   at nodes where displacement is constrained)

%--------------------------------------------------------------------------

%Error checks
if size(K,1) ~= size(K, 2)
    error('K must be square matrix');
end
if size(applied_forces, 2) ~= 2 | size(applied_displacements, 2) ~= 2
    error('applied_forces and applied_displacements must be two-column matrices');
end
if size(K,1) ~= numel(applied_forces) | size(K,1) ~= numel(applied_displacements)
    error('Number of elements in applied_forces and applied_displacements matrices must equal size of K');
end

%reshape applied_displacements and applied_forces to be vectors
applied_displacements = reshape(applied_displacements', size(K, 1),[]);
applied_forces = reshape(applied_forces', size(K, 1),[]);

%identify row/col associated with constrained disp (ii) and free disp (jj)
ii = find(~isnan(applied_displacements));
jj = find(isnan(applied_displacements));

u = zeros(size(K,1), 1);
f = zeros(size(K,1), 1);
u(jj) = K(jj, jj) \ (applied_forces(jj) - K(jj, ii) * applied_displacements(ii));
f(ii) = K(ii, jj) * u(jj) + K(ii, ii) * applied_displacements(ii);
u(ii) = applied_displacements(ii);

%reshape back to matrices
u = reshape(u, 2, [])';
f = reshape(f, 2, [])';
end