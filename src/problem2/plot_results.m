function plot_results(x_a, elem, u, Es, Ss, properties, result_dir, flag)
    % Plot and save FEA results for Mindlin Plate
    % Inputs:
    % x_a: Node coordinates
    % elem: Element connectivity
    % u: Global displacement vector [w, theta_x, theta_y, ...]
    % Es: Strains [kx, ky, kxy, gxz, gyz] per element
    % Ss: Stresses/Moments [Mx, My, Mxy, Qx, Qy] per element
    % properties: [E, nu, t, k]
    % result_dir: Directory to save images
    % flag: Element type
    
    [num_nodes, ~] = size(x_a);
    t = properties(3);
    
    % --- 1. Extract Nodal Results ---
    w = zeros(num_nodes, 1);
    theta_x = zeros(num_nodes, 1);
    theta_y = zeros(num_nodes, 1);
    
    for i = 1:num_nodes
        w(i) = u(3*i-2);
        theta_x(i) = u(3*i-1);
        theta_y(i) = u(3*i);
    end
    
    % Rotation Magnitude
    theta_mag = sqrt(theta_x.^2 + theta_y.^2);
    
    % --- 2. Extract Elemental Results (for flat shading) ---
    % Es: [kx, ky, kxy, gxz, gyz]
    % Ss: [Mx, My, Mxy, Qx, Qy]
    
    % Calculate Von Mises Stress at Top Surface (z = t/2)
    % sigma_x = 12*Mx*z/t^3 -> 6*Mx/t^2
    % sigma_y = 6*My/t^2
    % tau_xy  = 6*Mxy/t^2
    % sigma_vm = sqrt(sx^2 + sy^2 - sx*sy + 3*txy^2)
    
    factor = 6 / t^2;
    Mx = Ss(:, 1);
    My = Ss(:, 2);
    Mxy = Ss(:, 3);
    
    sx = factor * Mx;
    sy = factor * My;
    txy = factor * Mxy;
    
    sigma_vm = sqrt(sx.^2 + sy.^2 - sx.*sy + 3*txy.^2);
    
    % Calculate Curvature Magnitude (as a proxy for "Strain")
    % k_mag = sqrt(kx^2 + ky^2 + kxy^2)
    kx = Es(:, 1);
    ky = Es(:, 2);
    kxy = Es(:, 3);
    strain_mag = sqrt(kx.^2 + ky.^2 + kxy.^2);
    
    % --- 3. Plotting Helper Function ---
    function plot_nodal_contour(data, title_str, filename)
        f = figure('Visible', 'off'); % Don't show window
        patch('Faces', elem, 'Vertices', [x_a(:,1), x_a(:,2), zeros(size(data))], ...
              'FaceColor', 'interp', 'CData', data, 'EdgeColor', 'none');
        title(title_str);
        xlabel('X (m)'); ylabel('Y (m)');
        colorbar;
        axis equal;
        axis([0 1 0 1]);
        view(2);
        saveas(f, fullfile(result_dir, filename));
        close(f);
    end

    function plot_elemental_contour(data, title_str, filename)
        f = figure('Visible', 'off');
        % For elemental data, CData length must match number of faces (elements)
        patch('Faces', elem, 'Vertices', [x_a(:,1), x_a(:,2), zeros(num_nodes,1)], ...
              'FaceColor', 'flat', 'CData', data, 'EdgeColor', 'none');
        title(title_str);
        xlabel('X (m)'); ylabel('Y (m)');
        colorbar;
        axis equal;
        axis([0 1 0 1]);
        view(2);
        saveas(f, fullfile(result_dir, filename));
        close(f);
    end

    % --- 4. Generate Plots ---
    
    % 1. Deflection w
    plot_nodal_contour(w, 'Transverse Deflection w (m)', '1_Deflection_w.png');
    
    % 2. Rotation Magnitude
    plot_nodal_contour(theta_mag, 'Rotation Magnitude (rad)', '2_Rotation_Magnitude.png');
    
    % 3. Von Mises Stress (Top Surface)
    plot_elemental_contour(sigma_vm, 'Von Mises Stress (Top Surface) (Pa)', '3_Stress_VonMises.png');
    
    % 4. Curvature Magnitude (Strain proxy)
    plot_elemental_contour(strain_mag, 'Curvature Magnitude (1/m)', '4_Strain_Curvature.png');
    
    % 5. Deformed Shape (3D)
    f = figure('Visible', 'off');
    scale = 0.1 / max(abs(w)); % Scale to 10% of domain size
    if max(abs(w)) == 0; scale = 1; end
    
    patch('Faces', elem, 'Vertices', [x_a(:,1), x_a(:,2), w*scale], ...
          'FaceColor', 'interp', 'CData', w, 'EdgeColor', 'k', 'FaceAlpha', 0.8);
    title(['Deformed Shape (Scale: ' num2str(scale, '%.2e') ')']);
    xlabel('X'); ylabel('Y'); zlabel('w');
    colorbar;
    view(3);
    axis equal;
    grid on;
    saveas(f, fullfile(result_dir, '5_Deformed_Shape_3D.png'));
    close(f);

end
