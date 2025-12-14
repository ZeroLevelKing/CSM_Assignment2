
% Mesh generation
% Updated for Assignment 2 Problem 2
% Supports variable mesh density via input arguments

function [x_a,elem]=generate_mesh(flag, nx, ny)
     
    % Geometry dimensions (Square plate 1m x 1m)
    L = 1.0; 
    H = 1.0; 
    
    % Default mesh density if not provided
    if nargin < 3
        % Assignment requires:
        % Triangles (flag=1): 50, 200, 5000 elements -> nx=ny=5, 10, 50
        % Quads (flag=2): 25, 100, 2500 elements -> nx=ny=5, 10, 50
        nx = 10; 
        ny = 10;
    end
    
    dx = L / nx;
    dy = H / ny;
    
    % Initialize the nodal position
    num_nodes = (nx + 1) * (ny + 1);
    x_a = zeros(num_nodes, 2);
    
    node_idx = 1;
    for j = 0:ny
        for i = 0:nx
            x_a(node_idx, 1) = i * dx;
            x_a(node_idx, 2) = j * dy;
            node_idx = node_idx + 1;
        end
    end
    
    % Construct the connectivity table
    if flag == 2 % Quadrilateral elements
        num_elems = nx * ny;
        elem = zeros(num_elems, 4);
        el_idx = 1;
        for j = 1:ny
            for i = 1:nx
                % Node indices for the current element (clockwise)
                % n1 -- n2
                % |     |
                % n4 -- n3
                
                n1 = (j-1)*(nx+1) + i;
                n2 = n1 + 1;
                n3 = n2 + (nx+1);
                n4 = n1 + (nx+1);
                
                elem(el_idx, :) = [n1, n2, n3, n4];
                el_idx = el_idx + 1;
            end
        end
        
    elseif flag == 1 % Triangular elements
        num_elems = nx * ny * 2;
        elem = zeros(num_elems, 3);
        el_idx = 1;
        for j = 1:ny
            for i = 1:nx
                % Split quad into two triangles
                % Quad nodes:
                n1 = (j-1)*(nx+1) + i;
                n2 = n1 + 1;
                n3 = n2 + (nx+1);
                n4 = n1 + (nx+1);
                
                % Triangle 1: n1-n2-n4 (Top-Left)
                elem(el_idx, :) = [n1, n2, n4];
                el_idx = el_idx + 1;
                
                % Triangle 2: n2-n3-n4 (Bottom-Right)
                elem(el_idx, :) = [n2, n3, n4];
                el_idx = el_idx + 1;
            end
        end
    end
    
end