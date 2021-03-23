function measurement(x,rN, numDiodes)
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

    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];

    rB_3 = cross(rB[:,1], rB[:,2])   # This is a bit sketchy
    rB_3 = [rB rB_3]'; 
     
    # dy/dθ
    # rB_hat_1 = [0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    # rB_hat_2 = [0 -rB[3,2] rB[2,2]; rB[3,2] 0 -rB[1,2]; -rB[2,2] rB[1,2] 0];
    # dθ = [rB_hat_1; rB_hat_2];
    dθ = c .* n
    dθ = dθ * rB_3;
     
    # dy/dβ
    dβ = zeros(1, 3);

    # dy/dc     
    dc = n * rB_3

    # dy/dα 
    ndα = [-cos.(ϵ).*sin.(α) cos.(ϵ).*cos.(α) zeros(size(α))];
    dα = c .* ndα
    dα = dα * rB_3 

    # dy/dϵ 
    ndϵ = [-sin.(ϵ).*cos.(α) sin.(ϵ).*sin.(α) cos.(ϵ)]; # NOT sure why the middle term isnt negative 
    dϵ = c .* ndϵ
    dϵ = dϵ * rB_3


    
    # C = [rB_hat_1 zeros(3,3); 
    #      rB_hat_2 zeros(3,3)];
    # C = [rB_hat_1 zeros(3,3) zeros(3,numDiodes) zeros(3, numDiodes); 
    #      rB_hat_2 zeros(3,3) zeros(3,numDiodes) zeros(3, numDiodes)];

    C = [dθ; dβ; dc; dα; dϵ]


    y = rB[:];

    return y, C

end
