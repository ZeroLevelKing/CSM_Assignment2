% This function assembles the global stiffness matrix for Mindlin Plate
% input:
% elem: connectivity table
% x_a: coordinates of all the nodes
% properties: material property vector [E, nu, t, k]
% flag: 1 for Tri3, 2 for Quad4
% output:
% K: the global stiffness matrix
function [K]=K_matrix(elem,x_a,properties,flag)

    [num_nodes, ~] = size(x_a);
    [num_elem, num_nodes_per_elem] = size(elem);
    
    % Initialize global stiffness matrix
    % Use sparse for efficiency with large meshes
    K = sparse(3*num_nodes, 3*num_nodes);
    
    % Material properties
    E  = properties(1);
    nu = properties(2);
    t  = properties(3);
    k  = properties(4); % Shear correction factor
    
    % Constitutive Matrices
    % Bending Db
    % Db = E*t^3 / (12*(1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]
    Db_val = (E * t^3) / (12 * (1 - nu^2));
    Db = Db_val * [1,  nu, 0;
                   nu, 1,  0;
                   0,  0,  (1-nu)/2];
                   
    % Shear Ds
    % Ds = k*G*t * [1, 0; 0, 1]
    G = E / (2 * (1 + nu));
    Ds_val = k * G * t;
    Ds = Ds_val * [1, 0;
                   0, 1];
                   
    % Gauss Quadrature Rules
    if flag == 2 % Quad4
        % Full Integration (2x2) for Bending
        % xi, eta, weight
        val = 1/sqrt(3);
        g_pts_full = [-val, -val, 1.0;
                       val, -val, 1.0;
                       val,  val, 1.0;
                      -val,  val, 1.0];
                      
        % Reduced Integration (1x1) for Shear
        g_pts_red = [0.0, 0.0, 4.0];
        
    elseif flag == 1 % Tri3
        % 1-point integration (Centroid)
        % xi = 1/3, eta = 1/3, w = 0.5 (Area of master triangle is 0.5)
        g_pts_full = [1/3, 1/3, 0.5];
        g_pts_red  = [1/3, 1/3, 0.5]; 
    end
    
    % Loop over elements
    for e = 1:num_elem
        nodes = elem(e, :);
        xe = x_a(nodes, 1);
        ye = x_a(nodes, 2);
        
        % Element Stiffness Matrices
        Kb_e = zeros(3*num_nodes_per_elem, 3*num_nodes_per_elem);
        Ks_e = zeros(3*num_nodes_per_elem, 3*num_nodes_per_elem);
        
        % --- Bending Stiffness (Full Integration) ---
        for i = 1:size(g_pts_full, 1)
            xi = g_pts_full(i, 1);
            eta = g_pts_full(i, 2);
            w = g_pts_full(i, 3);
            
            [Bb, ~, detJ, ~] = B_matrix(xi, eta, xe, ye, flag);
            
            % K = integral(B' * D * B) dA
            % dA = detJ * dxi * deta (weight includes dxi*deta factors)
            Kb_e = Kb_e + Bb' * Db * Bb * detJ * w;
        end
        
        % --- Shear Stiffness (Reduced Integration) ---
        for i = 1:size(g_pts_red, 1)
            xi = g_pts_red(i, 1);
            eta = g_pts_red(i, 2);
            w = g_pts_red(i, 3);
            
            [~, Bs, detJ, ~] = B_matrix(xi, eta, xe, ye, flag);
            
            Ks_e = Ks_e + Bs' * Ds * Bs * detJ * w;
        end
        
        Ke = Kb_e + Ks_e;
        
        % Assemble into global K
        % Map local DOFs to global DOFs
        global_dofs = zeros(1, 3*num_nodes_per_elem);
        for i = 1:num_nodes_per_elem
            node_idx = nodes(i);
            global_dofs(3*i-2) = 3*node_idx - 2; % w
            global_dofs(3*i-1) = 3*node_idx - 1; % theta_x
            global_dofs(3*i)   = 3*node_idx;     % theta_y
        end
        
        K(global_dofs, global_dofs) = K(global_dofs, global_dofs) + Ke;
    end
end
