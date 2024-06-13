function [J_FIM,inv_J_FIM] = calculateRecursion(D_11_cont,D_12_cont,D_21_cont,D_22_cont, NUM_STATE_VARS) %#codegen
% Calculate Recursion


    ARR_SIZE = size(D_11_cont,3);
    

    J_FIM = zeros(NUM_STATE_VARS, NUM_STATE_VARS, ARR_SIZE);
    
    inv_J_FIM =zeros(NUM_STATE_VARS, NUM_STATE_VARS, ARR_SIZE);

    % J0 
    % TODO: Implement J0 as function argument
    J = zeros(NUM_STATE_VARS, NUM_STATE_VARS);

    
    
    for i = 1:ARR_SIZE
    
    
    
            J = D_22_cont(:,:,i) - (D_21_cont(:,:,i)*pinv(J + D_11_cont(:,:,i))) * D_12_cont(:,:,i);
        
            J_FIM(:,:,i) = J;
            
            inv_J_FIM(:,:,i) = pinv(J);
    
    
    end

end