function measurement(x, rN, numDiodes)
    # (1) Get Sun Vector (Right now I just make one up) 


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

    rs = B_Q_N * [1, 0, 0]; # For now, I am pretending that the sun vector is constant (i.e., satellite is rotating but stationary)
    rs_hat = [0    -rs[3]   rs[2];   rs[3]   0    -rs[1];  -rs[2] rs[1]     0];

    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];
     
    # dy/dθ
    rB_hat_1 = [0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    rB_hat_2 = [0 -rB[3,2] rB[2,2]; rB[3,2] 0 -rB[1,2]; -rB[2,2] rB[1,2] 0];
    dθ = c .* n
    dθ = dθ * rs_hat;  # THIS ISN'T RIGHT, BUT THE DIMENSIONS WORK FOR TESTING (should be sun vector)???
     
    # dy/dβ
    dβ = zeros(numDiodes, 3);

    # dy/dc     
    dc = n * rs; 

    # dy/dα 
    ndα = [-cos.(ϵ).*sin.(α) cos.(ϵ).*cos.(α) zeros(size(α))];
    dα = c .* ndα
    dα = dα * rs;

    # dy/dϵ 
    ndϵ = [-sin.(ϵ).*cos.(α) sin.(ϵ).*sin.(α) cos.(ϵ)]; # NOT sure why the middle term isnt negative 
    dϵ = c .* ndϵ
    dϵ = dϵ * rs;


    
    # C = [rB_hat_1 zeros(3,3); 
    #      rB_hat_2 zeros(3,3)];
    Cr = [rB_hat_1 zeros(3,3) zeros(3,3*numDiodes); # Jacobian for the vector measurement portion, [6, 6+3i]
         rB_hat_2 zeros(3,3) zeros(3,3*numDiodes)];

    Ci = [dθ dβ Diagonal(dc) Diagonal(dα) Diagonal(dϵ)]  # Jacobian for current measurement portion, [i, 6+3i]
    # C = [dθ dβ dc dα dϵ]  # Jacobian, [i, 9]

    C = [Cr; Ci];

    I_meas = c .* (n * rs) .+ 0; # Measured current, no albedo      

    y = [rB[:]; I_meas[:]]; 


    return y, C

end
