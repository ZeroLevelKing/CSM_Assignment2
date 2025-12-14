% Enforce the displacement boundary conditions in this function
% input:
% F: the original force vector
% K: the original stiffness matrix
% boundary: boolean vector (1 for constrained, 0 for free)
% dis: prescribed displacement values
% output:
% F: modified force vector
% K: modified stiffness matrix
function [F,K]=Enforce_BC(F,K,boundary,dis)

    [num_dofs, ~] = size(K);
    
    % Direct Substitution Method
    % Iterate over all degrees of freedom
    for i = 1:num_dofs
        if boundary(i) == 1
            % This DOF is constrained to value dis(i)
            prescribed_val = dis(i);
            
            % 1. Modify Force Vector to account for non-zero prescribed displacements
            % Subtract the column of K corresponding to the constrained DOF
            % multiplied by the prescribed value from the Force vector.
            % This moves the known terms to the RHS.
            % Note: We only need to do this for equations that are NOT constrained themselves,
            % but doing it for all is fine because we will overwrite F(i) later.
            if prescribed_val ~= 0
                F = F - K(:, i) * prescribed_val;
            end
            
            % 2. Zero out the row and column in K
            K(i, :) = 0;
            K(:, i) = 0;
            
            % 3. Set diagonal to 1
            K(i, i) = 1;
            
            % 4. Set the RHS force vector entry to the prescribed value
            F(i) = prescribed_val;
        end
    end

end