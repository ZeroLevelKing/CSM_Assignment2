% Calculate the global force vector for Gaussian distributed load
% input:
% x_a: coordinates of all the nodes
% elem: connectivity table
% LoadParams: [C, sigma, xc, yc]
% flag: 1 for Tri3, 2 for Quad4
% output:
% F: global force vector
function [F]=F_vector(x_a, elem, LoadParams, flag)
    
    C = LoadParams(1);
    sigma = LoadParams(2);
    xc = LoadParams(3);
    yc = LoadParams(4);
    
    [num_nodes, ~] = size(x_a);
    [num_elem, num_nodes_per_elem] = size(elem);
    
    F = zeros(3*num_nodes, 1);
    
    % Gauss Quadrature for Load Integration
    if flag == 2 % Quad4
        val = 1/sqrt(3);
        g_pts = [-val, -val, 1.0;
                  val, -val, 1.0;
                  val,  val, 1.0;
                 -val,  val, 1.0];
    elseif flag == 1 % Tri3
        % 1-point integration
        g_pts = [1/3, 1/3, 0.5];
    end
    
    for e = 1:num_elem
        nodes = elem(e, :);
        xe = x_a(nodes, 1);
        ye = x_a(nodes, 2);
        
        Fe = zeros(3*num_nodes_per_elem, 1);
        
        for i = 1:size(g_pts, 1)
            xi = g_pts(i, 1);
            eta = g_pts(i, 2);
            w = g_pts(i, 3);
            
            % Get shape functions and detJ
            [~, ~, detJ, N] = B_matrix(xi, eta, xe, ye, flag);
            
            % Calculate global coordinates of integration point
            x_g = N * xe;
            y_g = N * ye;
            
            % Evaluate Load
            R2 = (x_g - xc)^2 + (y_g - yc)^2;
            q = C * exp(-R2 / (2*sigma^2));
            
            % Integrate: Fe += N^T * q * detJ * w
            % We only add to w-DOFs (indices 1, 4, 7, ...)
            
            for n = 1:num_nodes_per_elem
                idx_w = 3*(n-1) + 1;
                Fe(idx_w) = Fe(idx_w) + N(n) * q * detJ * w;
            end
        end
        
        % Assemble
        for n = 1:num_nodes_per_elem
            node_idx = nodes(n);
            % Global indices
            idx_w_global = 3*node_idx - 2;
            
            % Local indices
            idx_w_local = 3*(n-1) + 1;
            
            F(idx_w_global) = F(idx_w_global) + Fe(idx_w_local);
        end
    end
end
