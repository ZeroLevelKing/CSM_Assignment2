% Set up the boundary conditions for Mindlin Plate (Assignment 2 Problem 2)
% Pinned edges on all 4 sides: w = 0
% input:
% x_a: the coordinates of all the nodes
% elem: connectivity table (not used for simple pinned BCs but kept for interface consistency)
% output:
% boundary: boolean flag for each DOF (3 per node), 1 for constrained, 0 for free
% dis: prescribed nodal displacement values
function [boundary,dis]=Boundary_conditions(x_a,elem)

    [num_nodes, ~] = size(x_a);
    
    % 3 DOFs per node: w, theta_x, theta_y
    boundary = zeros(3*num_nodes, 1);
    dis = zeros(3*num_nodes, 1);
    
    % Domain dimensions
    L = 1.0;
    H = 1.0;
    tol = 1e-5;
    
    % Loop over all nodes to check if they are on the boundary
    for i = 1:num_nodes
        x = x_a(i, 1);
        y = x_a(i, 2);
        
        % Check if node is on any of the 4 edges
        on_left   = abs(x - 0) < tol;
        on_right  = abs(x - L) < tol;
        on_bottom = abs(y - 0) < tol;
        on_top    = abs(y - H) < tol;
        
        if on_left || on_right || on_bottom || on_top
            % Pinned edge: Constrain w (DOF 1)
            % Global DOF index for w at node i is 3*i - 2
            idx_w = 3*i - 2;
            
            boundary(idx_w) = 1;
            dis(idx_w) = 0;
            
            % Note: For "Hard Support", rotations are free.
            % If "Soft Support" was required, we might constrain tangential rotation.
            % Instructions specify: "Assume Hard Support (w=0 only)"
        end
    end
end