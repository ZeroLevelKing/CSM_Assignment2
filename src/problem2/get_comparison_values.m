% Script to extract comparison values for report
clear; clc;

% Path to the result file
% Use absolute path based on script location to avoid CWD issues
current_file_path = mfilename('fullpath');
[current_dir, ~, ~] = fileparts(current_file_path);
result_file = fullfile(current_dir, '..', '..', 'result', 'problem2', 'Quad4_50x50', 'DATA.mat');

if ~exist(result_file, 'file')
    error('Result file not found at: %s\nPlease run main.m first to generate results.', result_file);
end

load(result_file);

% 1. Maximum Deflection (w)
% u structure: [w1, tx1, ty1, w2, tx2, ty2, ...]
w = u(1:3:end);
[max_w_mag, max_w_idx] = max(abs(w));
% The actual value might be negative (downward), so let's get the signed value
max_w_val = w(max_w_idx);

fprintf('--- Comparison Values (Quad4 50x50) ---\n');
fprintf('Maximum Deflection (w_max): %.6e m\n', max_w_val);

% 2. Maximum Von Mises Stress
% Ss is [num_elem x 5] -> [Mx, My, Mxy, Qx, Qy]
% Need thickness t from properties
% properties = [E, nu, t, k]
t = properties(3);

Mx = Ss(:, 1);
My = Ss(:, 2);
Mxy = Ss(:, 3);

% Calculate stress at top surface z = t/2
% sigma = M * z / I = M * (t/2) / (t^3/12) = 6*M/t^2
factor = 6 / t^2;

sx = factor * Mx;
sy = factor * My;
txy = factor * Mxy;

sigma_vm_elem = sqrt(sx.^2 + sy.^2 - sx.*sy + 3*txy.^2);

[max_vm_elem_val, max_vm_idx] = max(sigma_vm_elem);

fprintf('Maximum Von Mises Stress (Element Centroid):   %.6e Pa\n', max_vm_elem_val);

% --- Nodal Stress Smoothing (Better Estimate for Peak) ---
% Find the center node (0.5, 0.5)
center_node_idx = find(abs(x_a(:,1) - 0.5) < 1e-5 & abs(x_a(:,2) - 0.5) < 1e-5);

if ~isempty(center_node_idx)
    % Find elements sharing this node
    % elem is [num_elem x 4]
    [r, c] = find(elem == center_node_idx);
    connected_elems = r;
    
    if ~isempty(connected_elems)
        % Average the Von Mises stress from connected element centroids
        % This is a simple smoothing technique for constant/bilinear elements
        avg_vm = mean(sigma_vm_elem(connected_elems));
        fprintf('Maximum Von Mises Stress (Center Node Avg):  %.6e Pa\n', avg_vm);
    else
        fprintf('Warning: No elements found connected to center node.\n');
    end
else
    fprintf('Warning: Center node not found exactly.\n');
end

