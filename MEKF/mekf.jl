function mekf(x0, P0, W, V, rN, whist, yhist, dt)

    xhist = zeros(7,size(yhist,2)); # x = [q beta] = [q0 qi qj qk bx by bz]
    xhist[:,1] = x0;
    
    Phist = zeros(6,6,size(yhist,2));
    Phist[:,:,1] = P0;
    
    for k = 1:(size(yhist,2)-1)
       
        # Predict x, P
        x_p, A = prediction(xhist[:,k],whist[:,k],dt); # State prediction
        P_p = A*Phist[:,:,k]*A' + W; # Covariance prediction 
        
        yp, C = measurement(x_p[1:4],rN);         
        
        # Innovation
    
        z = yhist[:,k] - yp; 
        S = C*P_p*C' + V;
    
        # Kalman Gain
        L = P_p * C' * S^(-1); # ... / S;
        
        # Update
        dx = L*z; 
        dPhi = dx[1:3]; 
        dB = dx[4:6]; # Change in bias terms 
        
        theta_temp = (norm(dPhi));
        rTemp = dPhi / theta_temp; 
        
    #     dq = [0.5*dPhi; (1-0.125 * (dPhi'*dPhi))];
    
        # dq = [cos(theta_temp/2); rTemp*sin(theta_temp/2)];   ########################
        dq = [rTemp*sin(theta_temp/2); cos(theta_temp/2)];
        
        xhist[1:4,k+1] = qmult(x_p[1:4], dq); 
        xhist[5:7,k+1] = x_p[5:7] + dB; #bias update
        
        Phist[:,:,k+1] = (I(6) - L*C) * P_p * (I(6) - L*C)' + L*V*L';  #covariance update

    end
    
    return xhist, Phist
    
end
    