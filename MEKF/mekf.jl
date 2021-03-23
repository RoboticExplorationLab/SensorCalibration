function mekf(x0, P0, W, V, rN, whist, yhist, dt, numDiodes)

    xhist = zeros(size(x0,1),size(yhist,2)); # x = [q beta c α ϵ] = [q0 qi qj qk bx by bz] x [~1500]
    xhist[:,1] = x0;
    
    Phist = zeros(size(P0,1),size(P0,1),size(yhist,2));
    Phist[:,:,1] = P0;
    
    for k = 1:(size(yhist,2)-1)
       
        # Predict x, P
        x_p, A = prediction(xhist[:,k],whist[:,k],dt, numDiodes); # State prediction
        P_p = A*Phist[:,:,k]*A' + W; # Covariance prediction 
        
        yp, C = measurement(x_p,rN, numDiodes);         
        
        # Innovation
    
        z = yhist[:,k] - yp; 
        S = C*P_p*C' + V;   # Are the dimensions of C correct?
    
        # Kalman Gain
        L = P_p * C' * S^(-1); 
        
        # Update
        dx = L*z;   ###########
        dPhi = dx[1:3]; 
        dB = dx[4:6]; # Change in bias terms 
        dα = dx[7:(7+numDiodes)];
        dϵ = dx[(8+numDiodes):end];
        
        theta_temp = (norm(dPhi));
        rTemp = dPhi / theta_temp; 
        
        #     dq = [0.5*dPhi; (1-0.125 * (dPhi'*dPhi))];
    
        # dq = [cos(theta_temp/2); rTemp*sin(theta_temp/2)];   ########################
        dq = [rTemp*sin(theta_temp/2); cos(theta_temp/2)];

        # Why do we add these to x_p and not xhist? 
        xhist[1:4,k+1] = qmult(x_p[1:4], dq); 
        xhist[5:7,k+1] = x_p[5:7] + dB; #bias update
        xhist[8:(8+numDiodes),k+1] = x_p[8:(8+numDiodes)] + dα;
        xhist[(9+numDiodes):end,k+1] = x_p[(9+numDiodes):end] + dϵ;
        
        Phist[:,:,k+1] = (I(size(Phist,1)) - L*C) * P_p * (I(size(Phist,1)) - L*C)' + L*V*L';  #covariance update

    end
    
    return xhist, Phist
    
end
    