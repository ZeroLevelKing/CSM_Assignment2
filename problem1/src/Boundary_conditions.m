% Set up the boundary conditions
% input:
% x_a: the coordinates of all the nodes
% elem: connectivity table
% output:
% boundary: boolean flag for each nodal displacement in x and y direction, 1 for constrained, 0 for free
% dis: prescribed nodal displacement in x and y direction
% l_area: surface area for each node
function [boundary,dis,l_area]=Boundary_conditions(x_a,elem)

    
    % Displacements imposition
    % Note: this method only works for straight line boundaries
    % parallel to the coordinate axis
    % 1st col: axis label 1 for x, 2 for y
    % 2nd col: location
    % 3rd col: 1 for displacement in x, 2 for y, and 3 for in both x and y
    %          direction
    % 4th col: value of the prescribed displacement
    % 
    % e.g apply displacement boundary condition on the edge of the domain at
    % x=1.5, u(1.5, y)=(0.0 10.0)  ->  (1, 1.5, 2, 10.0)
        
        A(1,:)=[1 0 3 0];        
        [boundary,dis]=displa(x_a,A);
        
    % Forces imposition (area of the nodes) 
    % Note: this method only works for straight line boundaries
    % parallel to the coordinate axis
    % 1st col: axis label 1 for x, 2 for y
    % 2nd col: location
    % 
    % e.g apply traction boundary conditions to the edge of the domain at
    % x=0  ->  (1,0)
		%_ Traction is applied on the top horizontal, y=y_max
        B(1,:)=[2 max(x_a(:,2))];   
        [l_area]=dist(x_a,elem,B);
end

% Set up the displacement boundary conditions in this function
% input:
% x_a: coordinates of all the nodes
% A: prescribed displacement boundary condition
% output:
% boundary: boolean flag for each nodal displacement in x and y direction, 1 for constrained, 0 for free
% disp: prescribed nodal displacement in x and y direction
function [boundary,dis]=displa(x_a,A)
    
    [num_nodes, ~] = size(x_a);
    boundary = zeros(2*num_nodes, 1);
    dis = zeros(2*num_nodes, 1);
    
    tol = 1e-6;
    [num_bc, ~] = size(A);
    
    for k = 1:num_bc
        axis_idx = A(k, 1); % 1 for x, 2 for y
        loc = A(k, 2);
        dof_dir = A(k, 3); % 1=x, 2=y, 3=both
        val = A(k, 4);
        
        % Find nodes on the boundary line
        nodes_on_boundary = find(abs(x_a(:, axis_idx) - loc) < tol);
        
        for i = 1:length(nodes_on_boundary)
            node = nodes_on_boundary(i);
            
            % X-direction constraint
            if dof_dir == 1 || dof_dir == 3
                idx = 2*node - 1;
                boundary(idx) = 1;
                dis(idx) = val;
            end
            
            % Y-direction constraint
            if dof_dir == 2 || dof_dir == 3
                idx = 2*node;
                boundary(idx) = 1;
                dis(idx) = val;
            end
        end
    end

end

% Calculate the surface area associated with each node
% If the node is not a surface node and does not belong to the Neumann 
% boundary conditions, its surface area is initialized as 0
% input: 
% x_a: coordinates of all the nodes
% elem: connectivity table
% B: prescribed traction boundary condition
% output:
% l_area: the surface area associated to each node
function [l_area]=dist(x_a,elem,B)

    [num_nodes, ~] = size(x_a);
    l_area = zeros(num_nodes, 1);
    [num_elem, num_nodes_per_elem] = size(elem);
    [num_bc, ~] = size(B);
    tol = 1e-6;

    for k = 1:num_bc
        axis_idx = B(k, 1);
        loc = B(k, 2);
        
        for e = 1:num_elem
            nodes = elem(e, :);
            
            % Define edges based on element type
            if num_nodes_per_elem == 3
                edges = [1 2; 2 3; 3 1];
            else % Quad
                edges = [1 2; 2 3; 3 4; 4 1];
            end
            
            for j = 1:size(edges, 1)
                n1 = nodes(edges(j, 1));
                n2 = nodes(edges(j, 2));
                
                % Check if both nodes are on the boundary line
                if abs(x_a(n1, axis_idx) - loc) < tol && ...
                   abs(x_a(n2, axis_idx) - loc) < tol
               
                    % Calculate edge length
                    dx = x_a(n1, 1) - x_a(n2, 1);
                    dy = x_a(n1, 2) - x_a(n2, 2);
                    len = sqrt(dx^2 + dy^2);
                    
                    % Distribute half length to each node (lumped area)
                    l_area(n1) = l_area(n1) + len / 2;
                    l_area(n2) = l_area(n2) + len / 2;
                end
            end
        end
    end

end