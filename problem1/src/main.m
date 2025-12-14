

% main routine
% FEM code for 2D linear elasticity
% Li, B. 2010 April

clear

flag=1;   % 1: triangular element;   2: quadrilateral element.

% 1. Generatemesh, initialize the coordinates of nodes and element
%    connectivity table

    % x_a is a [Nx2] matrix, N is the total number of nodes
    % elem is a [Mx3] matrix for triangular elements and [Mx4] matrix for
    % quadrilateral elements, M is the total number of elements
    [x_a,elem]=generate_mesh(flag);
        
    % Check the initial mesh
    if flag==1        
        triplot(elem,x_a(:,1),x_a(:,2));
    elseif flag==2
        quadplot(elem,x_a(:,1),x_a(:,2));
    end
    
    % Area and geometric center
    % xg is a [Mx2] matrix, M is the total number of elements
    % Area is a vector with M entries
    [xg,Area]=g_center(x_a,elem);

    [nodes,dim]=size(x_a);
    [elements,NNE]=size(elem);
    
    
 % 2. Boundary conditions
    % auxiliaire data structure defined for the enforcement of boundary conditions
    % boundary: a vector with 2N entries, each element of the vector is a boolean flag 
    %           for the nodal displacement in x or y direction, 
    %           1 for constrained, 0 for free
    %           e.g. (1 1 0 0 1 0 ...) means node 1 has displacement boundary conditions
    %           applied in x and y direction, but node 3 has x displacement b.c.
    % disp    : a vector with 2N entries, each element of the vector is the prescribed 
    %           nodal displacement in x or y direction. For those free nodes, use (0 0)
    % l_area  : a vector with N entries, each element of the vector is the surface area
    %           associated with a node. If the node is not subjected to traction boundary
    %           conditions, its surface area is zero.
    [boundary,dis,l_area]=Boundary_conditions(x_a,elem);

% 3. Material Properties   
    E  = 2.1e9; % Young's modulus [Pa]
    nu = 0.3;   % Poisson ratio

    properties(1)=E;
    properties(2)=nu;

% 4. B matrix: B is a cell structure where B(i) is the B matrix of element i    
    [B, N]=B_matrix(x_a,elem,xg,Area,flag);

% 5. K matrix: K is a 2Nx2N matrix, N is the total number of nodes
    [K]=K_matrix(B,elem,x_a,Area,properties);

% 6. Forces vector: F is a vector with 2N entries
    Load = [0, -2e4]; % traction [N/m] 
    [F]=F_vector(x_a,Load,l_area);

% 7. Enforce Essential Boundary Condition
    [F,K]=Enforce_BC(F,K,boundary,dis,x_a);
    
% 8. Solve the problem
    [u]=K\F;
    
% 9. Compute the strain, stress and pressure from the displacement
    [Es,Ss,P]=constitutive(B,properties,u,elem,dim);

% 10. Post-processor    
    h=10; % amplification factor
    
    % Create result directory if it doesn't exist
    result_dir = fullfile(pwd, 'result');
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end

    plot_results(x_a,elem,P,Ss,u,h,N,flag,result_dir);

    save(fullfile(result_dir, 'DATA.mat'), 'Es', 'Ss', 'P', 'u');














