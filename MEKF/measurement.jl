function measurement(x, rN, numDiodes)
    # Generates the "measured" body vectors for a set of 
    #      newtonian (inertial) vectors using quaternion q
    #      (What our measurement would be given our estimated attitude)
    # Inputs:
    #   - q: quaternion representing current orientation.
    #           Scalar first (10, \vec{q})                             |  [4 x 1]
    #   - rN: Pair of unit vectors in the newtonian (inertial) frame   |  [3 x 2]
    # Outputs:
    #   - y: Pair of unit vectors in body frame corresponding to rN & q
    #           (Flattened)                                            |  [6 x 1]
    #   - C: ()
    #           dy/dx = [dy/dPhi; dy/dBeta] = [dy/dPhi; 0]             |  [6 x 6]

    x = x[:]
    q = x[1:4]
    β = x[5:7]
    c = x[8:(7+numDiodes)]  # Find better way to do these
    α = x[(8+numDiodes):(7+2*numDiodes)]
    ϵ = x[(8+2*numDiodes):end]

    # DCM from quaternion (flipped)
    B_Q_N = dcm_from_q(q)';                           # # # #  # # # #  # # # #  

    # this is what the measurement would be given our estimated attitude
    rB = B_Q_N*rN;

    rs = B_Q_N * rN[:,1]; # Assume that rN[:,1] is the sun vector %%%%%%%%%%%%%%
    rs_hat = [0     -rs[3]   rs[2];   
              rs[3]    0    -rs[1];  
             -rs[2]  rs[1]     0];

    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];

    # theta = n[1,:]' * rs[:,1];
    # theta = theta / (norm(n[1,:])*norm(rs[:,1]))
    # theta = acos(theta) * 180 / pi

    # if theta > 90
    #     println("Error - greater than 90")
    # elseif theta < 0
    #     println("Error - theta < 0")
    # end

     
    # # dy/dθ
    rB_hat_1 = [0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    rB_hat_2 = [0 -rB[3,2] rB[2,2]; rB[3,2] 0 -rB[1,2]; -rB[2,2] rB[1,2] 0];
    dθ = c .* n
    dθ = dθ * rs_hat;  
     
    # dy/dβ
    dβ = zeros(numDiodes, 3);

    # dy/dc     
    dc = n * rs; 

    # dy/dα 
    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    dα = c .* (ndα * rs);

    # dy/dϵ 
    ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term)
    # ndϵ = [(-sin.(ϵ).*cos.(α)) (sin.(ϵ).*sin.(α)) cos.(ϵ)]; # NOT sure why the middle term isnt negative 
    dϵ = c .* (ndϵ * rs);

    
    # C = [rB_hat_1 zeros(3,3); 
    #      rB_hat_2 zeros(3,3)];
    Cr = [rB_hat_1 zeros(3,3) zeros(3,3*numDiodes); # Jacobian for the vector measurement portion, [6, 6+3i]
         rB_hat_2 zeros(3,3) zeros(3,3*numDiodes)];

    Ci = [dθ dβ Diagonal(dc) Diagonal(dα) Diagonal(dϵ)]  # Jacobian for current measurement portion, [i, 6+3i]
    
    
    # Cr = [rB_hat_1 zeros(3,3) zeros(3,3); # Jacobian for the vector measurement portion, [6, 9]
    #       rB_hat_2 zeros(3,3) zeros(3,3)];
    # Ci = [dθ dβ dc dα dϵ]  # Jacobian, [i, 9]

    C = [Cr; Ci];

    I_meas = c .* (n * rs) .+ 0; # Measured current, no albedo      

    y = [rB[:]; I_meas[:]]; 

    # y = rB[:];

    return y, C

end
