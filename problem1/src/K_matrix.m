% This function assembles the global stiffness matrix
% input:
% B: the B matrix for all the elements
% elem: connectivity table
% x_a: coordinates of all the nodes
% jacobians: jacobians of all the elements
% properties: material property vector
% output:
% K: the global stiffness matrix
function [K]=K_matrix(B,elem,x_a,jacobians,properties)

    [num_nodes, ~] = size(x_a);
    [num_elem, num_nodes_per_elem] = size(elem);
    
    % Initialize global stiffness matrix
    K = zeros(2*num_nodes, 2*num_nodes);
    
    % Material properties
    E = properties(1);
    nu = properties(2);
    
    % Plane Stress D matrix
    % D = E / (1-nu^2) * [1  nu 0;
    %                     nu 1  0;
    %                     0  0  (1-nu)/2];
    factor = E / (1 - nu^2);
    D = factor * [1,  nu, 0;
                  nu, 1,  0;
                  0,  0,  (1-nu)/2];
              
    % Loop over elements
    for e = 1:num_elem
        Be = B{e};          % B matrix for element e [3 x 2*NNE]
        Ae = jacobians(e);  % Area of element e (passed as 'jacobians' argument)
        
        % Element stiffness matrix
        % Ke = integral(B' * D * B) dV
        % Assuming constant B and unit thickness: Ke = B' * D * B * Area
        Ke = Be' * D * Be * Ae;
        
        % Assemble into global K
        nodes = elem(e, :);
        
        % Map local DOFs to global DOFs
        % Local:  1, 2, 3, 4, ...
        % Global: 2*n1-1, 2*n1, 2*n2-1, 2*n2, ...
        
        global_dofs = zeros(1, 2*num_nodes_per_elem);
        for i = 1:num_nodes_per_elem
            node_idx = nodes(i);
            global_dofs(2*i-1) = 2*node_idx - 1; % x-dof
            global_dofs(2*i)   = 2*node_idx;     % y-dof
        end
        
        % Add Ke to K
        K(global_dofs, global_dofs) = K(global_dofs, global_dofs) + Ke;
    end
  
end

