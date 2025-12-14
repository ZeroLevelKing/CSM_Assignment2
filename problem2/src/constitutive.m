% Calculate strains and stresses (Moments/Shears) at element centroids
% input:
% elem: connectivity table
% x_a: coordinates of all the nodes
% properties: material properties
% u: global displacement vector
% flag: element type
% output:
% Es: Strains [kx, ky, kxy, gxz, gyz] per element
% Ss: Stresses [Mx, My, Mxy, Qx, Qy] per element
function [Es,Ss]=constitutive(elem, x_a, properties, u, flag)
    
    [num_elem, num_nodes_per_elem] = size(elem);
    
    % Material properties
    E  = properties(1);
    nu = properties(2);
    t  = properties(3);
    k  = properties(4);
    
    % Constitutive Matrices
    Db_val = (E * t^3) / (12 * (1 - nu^2));
    Db = Db_val * [1,  nu, 0;
                   nu, 1,  0;
                   0,  0,  (1-nu)/2];
                   
    G = E / (2 * (1 + nu));
    Ds_val = k * G * t;
    Ds = Ds_val * [1, 0;
                   0, 1];
                   
    % Centroid integration point
    if flag == 2 % Quad4
        xi = 0; eta = 0;
    elseif flag == 1 % Tri3
        xi = 1/3; eta = 1/3;
    end
    
    % Initialize outputs
    % Es: [kx; ky; kxy; gxz; gyz] for each element (5 components)
    % Ss: [Mx; My; Mxy; Qx; Qy] for each element (5 components)
    Es = zeros(num_elem, 5);
    Ss = zeros(num_elem, 5);
    
    for e = 1:num_elem
        nodes = elem(e, :);
        xe = x_a(nodes, 1);
        ye = x_a(nodes, 2);
        
        % Get element displacements
        ue = zeros(3*num_nodes_per_elem, 1);
        for i = 1:num_nodes_per_elem
            node_idx = nodes(i);
            ue(3*i-2) = u(3*node_idx-2);
            ue(3*i-1) = u(3*node_idx-1);
            ue(3*i)   = u(3*node_idx);
        end
        
        % Get B matrices at centroid
        [Bb, Bs, ~, ~] = B_matrix(xi, eta, xe, ye, flag);
        
        % Calculate Curvatures and Shear Strains
        kappa = Bb * ue; % [kx; ky; kxy]
        gamma = Bs * ue; % [gxz; gyz]
        
        % Calculate Moments and Shear Forces
        Moments = Db * kappa; % [Mx; My; Mxy]
        Shears  = Ds * gamma; % [Qx; Qy]
        
        Es(e, :) = [kappa', gamma'];
        Ss(e, :) = [Moments', Shears'];
    end
end
