function mekf(x0, P0, W, V, rN, whist, yhist, dt, numDiodes, eclipse)

    xhist = zeros(size(x0,1),size(yhist,2)); # x = [q β c α ϵ] = [qi qj qk q0  bx by bz  ci  αi  ϵi ]   | [7 + 3i x n]
    xhist[:,1] = x0;
    
    Phist = zeros(size(P0,1),size(P0,1),size(yhist,2));
    Phist[:,:,1] = P0;
    
    for k = 1:(size(yhist,2)-1)
       
        # Predict x, P
        x_p, A = prediction(xhist[:,k],whist[:,k],dt, numDiodes); # State prediction
        P_p = A*Phist[:,:,k]*A' + W; # Covariance prediction 

        # Measurement
        z = [] 
        C = Array{Float64}(undef, 0, size(P_p,1))
        Vs = Array{Float64}(undef)  ############ ASSUMES NO CORRELATION AT ALL!!!!

        if eclipse[k] > 0  #  Sun Vector measurement
            sN = rN[1:3,k];
            V_sun = V[1:3, 1:3]; 
            yp_sun, C_sun = sun_measurement(x_p, sN, numDiodes) #, eclipse[k])
            z_sun = yhist[1:3,k]  - yp_sun
            z = [z[:]; z_sun[:]]
            C = [C; C_sun]
            Vs = [Vs[:]; diag(V_sun)[:]]
        end

        if true        #  Magnetic Field measurement
            bN = rN[4:6,k];
            V_mag = V[4:6, 4:6];
            yp_mag, C_mag = mag_measurement(x_p, bN, numDiodes)
            z_mag = yhist[4:6,k] - yp_mag
            z = [z[:]; z_mag[:]]
            C = [C; C_mag];
            # V_mag = [V_mag[1,1], V_mag[2,2], V_mag[3,3]]
            Vs = [Vs[:]; diag(V_mag)[:]]
        end
        
        if eclipse[k] > 0   #  Diode Current Measurement
            V_cur = V[7:end, 7:end];
            yp_cur, C_cur = current_measurement(x_p, sN, numDiodes, eclipse[k])
            z_cur = yhist[7:end,k] - yp_cur       
            z = [z[:]; z_cur[:]]
            C = [C; C_cur]
            Vs = [Vs[:]; diag(V_cur)[:]]
        end



        # Innovation   
        Vs = Vs[2:end] # Get rid of the initial 0 term   
        Vk = Diagonal(Vs) 
        S = C*P_p*C' + Vk;   
    
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
        
        Phist[:,:,k+1] = (I(size(Phist,1)) - L*C) * P_p * (I(size(Phist,1)) - L*C)' + L*Vk*L';  
    end
    
    
    return xhist, Phist
    
end
