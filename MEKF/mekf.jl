function mekf(x0, P0, W, V, rN, whist, yhist, dt, numDiodes, eclipse)

    xhist = zeros(size(x0,1),size(yhist,2)); # x = [q beta c α ϵ] = [q0 qi qj qk bx by bz ci αi ϵi ] x [N]
    xhist[:,1] = x0;
    
    Phist = zeros(size(P0,1),size(P0,1),size(yhist,2));
    Phist[:,:,1] = P0;
    
    for k = 1:(size(yhist,2)-1)
       
        # Predict x, P
        x_p, A = prediction(xhist[:,k],whist[:,k],dt, numDiodes); # State prediction
        P_p = A*Phist[:,:,k]*A' + W; # Covariance prediction 
        
        
        rNVecs = [rN[1:3,k] rN[4:6,k]]; #######################
        # rNVecs = [rN[1:3] rN[4:6]] ##############################
        yp, C = measurement(x_p, rNVecs, numDiodes, eclipse[k]);   # Gives us what we would expect to get out of yhist given our predicted x
        
        # Innovation
        z = yhist[:,k] - yp[:];
        S = C*P_p*C' + V;   
    
        # Kalman Gain
        L = P_p * C' * S^(-1); 
        
        # Update
        dx = L*z;          
        dPhi = dx[1:3]; 
        drest = dx[4:end]

        theta_temp = (norm(dPhi));
        rTemp = dPhi / theta_temp; 
        
        dq = [rTemp*sin(theta_temp/2); cos(theta_temp/2)];

        xhist[1:4,k+1] = qmult(x_p[1:4], dq); 
        xhist[5:end, k+1] = x_p[5:end] + drest;
        
        Phist[:,:,k+1] = (I(size(Phist,1)) - L*C) * P_p * (I(size(Phist,1)) - L*C)' + L*V*L';  #covariance update
    end
    
    return xhist, Phist
    
end