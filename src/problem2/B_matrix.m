% Calculate the B matrices (Bending and Shear) at a specific integration point
%
% input:
% xi, eta : parametric coordinates of the integration point
% xe, ye  : coordinates of the element nodes (vectors)
% flag    : 1 for triangular elements (Tri3), 2 for quadrilateral elements (Quad4)
%
% output:
% Bb   : Bending B-matrix [3 x 3*num_nodes]
% Bs   : Shear B-matrix [2 x 3*num_nodes]
% detJ : Determinant of the Jacobian
% N    : Shape function values at (xi, eta)
function [Bb, Bs, detJ, N]=B_matrix(xi, eta, xe, ye, flag)

    if flag == 2 % Quad4
        % Node local coordinates (Standard isoparametric Quad4)
        % 4 -- 3
        % |    |
        % 1 -- 2
        xi_nodes  = [-1,  1,  1, -1];
        eta_nodes = [-1, -1,  1,  1];
        
        % Shape functions and derivatives wrt xi, eta
        N = zeros(1,4);
        dNdxi = zeros(1,4);
        dNdeta = zeros(1,4);
        
        for i = 1:4
            N(i) = 0.25 * (1 + xi_nodes(i)*xi) * (1 + eta_nodes(i)*eta);
            dNdxi(i) = 0.25 * xi_nodes(i) * (1 + eta_nodes(i)*eta);
            dNdeta(i) = 0.25 * eta_nodes(i) * (1 + xi_nodes(i)*xi);
        end
        
    elseif flag == 1 % Tri3
        % Standard triangle: (0,0), (1,0), (0,1)
        % N1 = 1 - xi - eta
        % N2 = xi
        % N3 = eta
        
        N = [1-xi-eta, xi, eta];
        dNdxi = [-1, 1, 0];
        dNdeta = [-1, 0, 1];
    end
    
    % Jacobian Matrix
    % J = [dx/dxi  dy/dxi;
    %      dx/deta dy/deta]
    J = [dNdxi * xe,  dNdxi * ye;
         dNdeta * xe, dNdeta * ye];
     
    detJ = det(J);
    
    % Inverse Jacobian
    % Handle singular Jacobian if necessary (mesh distortion)
    if abs(detJ) < 1e-10
        invJ = zeros(2,2); % Should trigger error or warning in robust code
    else
        invJ = inv(J);
    end
    
    % Derivatives wrt global coordinates x, y
    % [dN/dx] = invJ * [dN/dxi]
    % [dN/dy]          [dN/deta]
    
    dN_dxi_eta = [dNdxi; dNdeta];
    dN_dx_y = invJ * dN_dxi_eta;
    
    dNdx = dN_dx_y(1, :);
    dNdy = dN_dx_y(2, :);
    
    % Initialize B matrices
    num_nodes = length(xe);
    Bb = zeros(3, 3*num_nodes);
    Bs = zeros(2, 3*num_nodes);
    
    % Construct B matrices
    % DOFs per node: [w, theta_x, theta_y]
    % theta_x: Rotation about y-axis -> u = z*theta_x
    % theta_y: Rotation about x-axis -> v = z*theta_y
    
    for i = 1:num_nodes
        % Column indices for node i
        idx = 3*(i-1) + (1:3); 
        
        % Bending Strain (Curvatures)
        % kx  = d(theta_x)/dx
        % ky  = d(theta_y)/dy
        % kxy = d(theta_x)/dy + d(theta_y)/dx
        
        Bb(:, idx) = [0, dNdx(i), 0;
                      0, 0,       dNdy(i);
                      0, dNdy(i), dNdx(i)];
                      
        % Shear Strain
        % gxz = dw/dx + theta_x
        % gyz = dw/dy + theta_y
        
        Bs(:, idx) = [dNdx(i), N(i), 0;
                      dNdy(i), 0,    N(i)];
    end

end
