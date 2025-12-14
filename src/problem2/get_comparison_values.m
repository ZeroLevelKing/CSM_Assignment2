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

sigma_vm = sqrt(sx.^2 + sy.^2 - sx.*sy + 3*txy.^2);

[max_vm, max_vm_idx] = max(sigma_vm);

fprintf('Maximum Von Mises Stress:   %.6e Pa\n', max_vm);
