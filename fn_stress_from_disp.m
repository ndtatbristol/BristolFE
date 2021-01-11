function [sigma_11, sigma_22, sigma_12] = fn_stress_from_disp(Q, u)
%SUMMARY
%   Converts nodal displacements to stress within elements
%INPUTS
%   Q - m x n sparse matrix of conversion coefficients
%   u - n/2 x 2 matrix of nodal displacements
%OUTPUTS
%   sigma_11, sigma_22, sigma_12 - m/3 x 3 vectors of stress components in
%   the m elements

%--------------------------------------------------------------------------
%Error checks
if size(u, 2) ~= 2
    error('u must be two column matrix');
end
if size(Q, 2) ~= numel(u)
    error('Number of columns in Q must equal number of elements in u');
end

%reshape u
u = reshape(u', size(Q, 2),[]);
    
stress = Q * u(:);

%reshape stress to get components into columns
stress = reshape(stress, 3, [])';
% sigma = zeros(max(Q_matrix_el), 3);
% for ii = 1:3
%     jj = find(Q_matrix_stress_dir == ii);
%     sigma(Q_matrix_el(jj), ii) = stress(jj);
% end
sigma_11 = stress(:,1);
sigma_22 = stress(:,2);
sigma_12 = stress(:,3);
end