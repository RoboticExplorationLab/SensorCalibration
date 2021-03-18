function measurement(q,rN)
    # Generates the "measured" body vectors for a set of 
    #      newtonian (inertial) vectors using quaternion q
    #      (What our measurement would be given our estimated attitude)
    # Inputs:
    #   - q: quaternion representing current orientation.
    #           Scalar first (10, \vec{q})                              |  [4 x 1]
    #   - rN: Pair of unit vectors in the newtonian (inertial) frame   |  [3 x 2]
    # Outputs:
    #   - y: Pair of unit vectors in body frame corresponding to rN & q
    #           (Flattened)                                            |  [6 x 1]
    #   - C: ()
    #           dy/dx = [dy/dPhi; dy/dBeta] = [dy/dPhi; 0]             |  [6 x 6]

    # DCM from quaternion (flipped)
    B_Q_N = dcm_from_q(q)';                           # # # #  # # # #  # # # #  

    # this is what the measurement would be given our estimated attitude
    rB = B_Q_N*rN;
     
    rB_hat_1 = [0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    rB_hat_2 = [0 -rB[3,2] rB[2,2]; rB[3,2] 0 -rB[1,2]; -rB[2,2] rB[1,2] 0];
    
    C = [rB_hat_1 zeros(3,3); 
         rB_hat_2 zeros(3,3)];

    y = rB[:];

    return y, C

end
