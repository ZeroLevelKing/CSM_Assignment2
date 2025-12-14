

function plot_results(x_a,elem,P,Ss,d,h,N,flag,result_dir)

    [elements,NNE]=size(elem);
    [nodes,sp]=size(x_a);
    
    % Extract displacements
    ux = d(1:2:end);
    uy = d(2:2:end);
    u_mag = sqrt(ux.^2 + uy.^2);
    
    % Deformed coordinates
    x_def = x_a(:,1) + h * ux;
    y_def = x_a(:,2) + h * uy;
    
    % --- 1. Displacement Cloud Map ---
    figure('Name', 'Displacement Magnitude', 'Visible', 'off');
    patch('Faces', elem, 'Vertices', [x_def, y_def], 'FaceVertexCData', u_mag, ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
    colorbar;
    title(sprintf('Displacement Magnitude (Amplification: %d)', h));
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    colormap jet;
    
    % Save Displacement Plot
    saveas(gcf, fullfile(result_dir, 'displacement_contour.png'));
    saveas(gcf, fullfile(result_dir, 'displacement_contour.fig'));
    
    % --- 2. Stress Cloud Map (Von Mises) ---
    % Calculate Von Mises stress for each element
    vm_stress_elem = zeros(elements, 1);
    for e = 1:elements
        % Extract stress components for element e
        % Ss is stored as [sig_x, sig_y, tau_xy] for each element
        idx = (e-1)*3;
        sig_x = Ss(idx+1);
        sig_y = Ss(idx+2);
        tau_xy = Ss(idx+3);
        
        % Von Mises for Plane Stress
        vm_stress_elem(e) = sqrt(sig_x^2 - sig_x*sig_y + sig_y^2 + 3*tau_xy^2);
    end
    
    % Smooth element stresses to nodes
    vm_stress_nodes = zeros(nodes, 1);
    node_count = zeros(nodes, 1);
    
    for e = 1:elements
        for j = 1:NNE
            node_idx = elem(e, j);
            vm_stress_nodes(node_idx) = vm_stress_nodes(node_idx) + vm_stress_elem(e);
            node_count(node_idx) = node_count(node_idx) + 1;
        end
    end
    
    % Average
    vm_stress_nodes = vm_stress_nodes ./ max(node_count, 1);
    
    figure('Name', 'Von Mises Stress', 'Visible', 'off');
    patch('Faces', elem, 'Vertices', [x_def, y_def], 'FaceVertexCData', vm_stress_nodes, ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
    colorbar;
    title('Von Mises Stress (Pa)');
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    colormap jet;
    
    % Save Stress Plot
    saveas(gcf, fullfile(result_dir, 'stress_contour.png'));
    saveas(gcf, fullfile(result_dir, 'stress_contour.fig'));
    
    % --- 3. Original Deformed Mesh ---
    figure('Name', 'Deformed Mesh', 'Visible', 'off');
    patch('Faces', elem, 'Vertices', x_a, 'FaceColor', 'none', 'EdgeColor', 'b');
    hold on;
    patch('Faces', elem, 'Vertices', [x_def, y_def], 'FaceColor', 'none', 'EdgeColor', 'r');
    title('Deformed Mesh (Red) vs Original (Blue)');
    axis equal;
    
    saveas(gcf, fullfile(result_dir, 'deformed_mesh.png'));
    
    % Close all figures
    close all;

end