% main routine
% FEM code for Mindlin Plate Bending (Assignment 2, Problem 2)
% Refactored for 3 DOFs per node (w, theta_x, theta_y)
% Supports Mesh Convergence Study for both Tri3 and Quad4 elements

clear
close all; 

% Define Mesh Densities for Convergence Study
% Case 1: 5x5 (Coarse)
% Case 2: 10x10 (Medium)
% Case 3: 50x50 (Fine)
mesh_cases = [5, 10, 50];

% Element types to run: 1=Tri3, 2=Quad4
element_types = [1, 2];

% Material Properties   
E  = 210e9;   % Young's modulus [Pa]
nu = 0.3;     % Poisson ratio
t  = 0.05;    % Thickness [m]
k  = 5/6;     % Shear correction factor

properties(1)=E;
properties(2)=nu;
properties(3)=t;
properties(4)=k;

% Load Parameters
% Amplitude is negative (-1e6) to apply a downward force (in -z direction)
% Assignment Formula: F = 10^6 * (1/(sigma*sqrt(2*pi))) * exp(...)
sigma = 0.1;
Amp = -1e6 * (1 / (sigma * sqrt(2*pi)));
LoadParams = [Amp, sigma, 0.5, 0.5]; 

for type_idx = 1:length(element_types)
    flag = element_types(type_idx);
    
    if flag == 1
        type_str = 'Tri3';
    else
        type_str = 'Quad4';
    end
    
    fprintf('================================================\n');
    fprintf('Starting Analysis for Element Type: %s\n', type_str);
    
    for i = 1:length(mesh_cases)
        nx = mesh_cases(i);
        ny = nx;
        
        fprintf('------------------------------------------------\n');
        fprintf('Running Case: %s, Mesh %dx%d\n', type_str, nx, ny);
        
        % 1. Generate Mesh
        [x_a, elem] = generate_mesh(flag, nx, ny);
        [nodes, ~] = size(x_a);
        [elements, ~] = size(elem);
        fprintf('Number of Nodes: %d, Number of Elements: %d\n', nodes, elements);
        
        % 2. Boundary Conditions
        [boundary, dis] = Boundary_conditions(x_a, elem);
        
        % 3. Stiffness Matrix
        fprintf('Assembling Stiffness Matrix...\n');
        [K] = K_matrix(elem, x_a, properties, flag);
        
        % 4. Force Vector
        fprintf('Assembling Force Vector...\n');
        [F] = F_vector(x_a, elem, LoadParams, flag);
        
        % 5. Enforce BCs
        [F, K] = Enforce_BC(F, K, boundary, dis);
        
        % 6. Solve
        fprintf('Solving System...\n');
        [u] = K \ F;
        
        % 7. Post-processing
        fprintf('Computing Stresses...\n');
        [Es, Ss] = constitutive(elem, x_a, properties, u, flag);
        
        % 8. Save Results
        case_name = sprintf('%s_%dx%d', type_str, nx, ny);
        % Save results to ../../result/problem2 relative to this script
        current_file_path = mfilename('fullpath');
        [current_dir, ~, ~] = fileparts(current_file_path);
        result_dir = fullfile(current_dir, '..', '..', 'result', 'problem2', case_name);
        
        if ~exist(result_dir, 'dir')
            mkdir(result_dir);
        end
        
        plot_results(x_a, elem, u, Es, Ss, properties, result_dir, flag);
        
        save_path = fullfile(result_dir, 'DATA.mat');
        save(save_path, 'Es', 'Ss', 'u', 'x_a', 'elem', 'properties', 'nx', 'ny', 'flag');
        
        fprintf('Results saved to: %s\n', result_dir);
    end
end

fprintf('================================================\n');
fprintf('All analysis cases completed.\n');
