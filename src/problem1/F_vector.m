% Assemble the force vector in this function
% input:
% x_a    : nodal vector
% Load   : traction on each node
% l_area : surface area associated with each node
% output:
% F      : the global force vector
function [F]=F_vector(x_a,Load,l_area)

    [num_nodes, ~] = size(x_a);
    F = zeros(2*num_nodes, 1);
    
    % Load is a vector [Tx, Ty] (traction in N/m)
    Tx = Load(1);
    Ty = Load(2);
    
    for i = 1:num_nodes
        % Only apply force if the node has a non-zero associated boundary length
        if l_area(i) > 0
            % Force = Traction * Length
            Fx = Tx * l_area(i);
            Fy = Ty * l_area(i);
            
            % Add to global force vector
            idx_x = 2*i - 1;
            idx_y = 2*i;
            
            F(idx_x) = Fx;
            F(idx_y) = Fy;
        end
    end
    
end



    