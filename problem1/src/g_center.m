% Calculate the barycenter and surface area of each element
%
% input:
% x_a  : coordinates of all the nodes
% elem : connectivity table
%
% output:
% xg   : barycenters of all the elements
% Area : surface areas of all the elements
function [xg,Area]=g_center(x_a,elem)

    [num_elem, num_nodes_per_elem] = size(elem);
    xg = zeros(num_elem, 2);
    Area = zeros(num_elem, 1);

    for i = 1:num_elem
        % Get coordinates of nodes for the current element
        nodes_idx = elem(i, :);
        x = x_a(nodes_idx, 1);
        y = x_a(nodes_idx, 2);
        
        % Calculate Barycenter (Geometric Center)
        xg(i, 1) = mean(x);
        xg(i, 2) = mean(y);
        
        % Calculate Area
        if num_nodes_per_elem == 3
            % Triangle
            % Area = 0.5 * |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)|
            Area(i) = 0.5 * abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)));
            
        elseif num_nodes_per_elem == 4
            % Quadrilateral
            % Using Shoelace formula (Surveyor's formula)
            % Area = 0.5 * |sum(xi*y(i+1) - x(i+1)*yi)|
            % Note: Indices wrap around (4 -> 1)
            
            term1 = x(1)*y(2) - x(2)*y(1);
            term2 = x(2)*y(3) - x(3)*y(2);
            term3 = x(3)*y(4) - x(4)*y(3);
            term4 = x(4)*y(1) - x(1)*y(4);
            
            Area(i) = 0.5 * abs(term1 + term2 + term3 + term4);
        end
    end
        
end