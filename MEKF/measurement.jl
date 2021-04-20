function g(x, rN, numDiodes)
    x = x[:]
    q = x[1:4]
    β = x[5:7]

    # DCM from quaternion (flipped)
    B_Q_N = dcm_from_q(q)';                           

    # this is what the measurement would be given our estimated attitude
    rB = B_Q_N*rN;

    rs = B_Q_N * rN[:,1]; # Assumes that rN[:,1] is the sun vector %%%%%%%%%%%%%%
    rs_hat = [0     -rs[3]   rs[2];   rs[3]    0    -rs[1];  -rs[2]  rs[1]     0];

    c = x[8:(7+numDiodes)]  # Find better way to do these
    α = x[(8+numDiodes):(7+2*numDiodes)]
    ϵ = x[(8+2*numDiodes):end]

    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];

    rB_hat_1 = [0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    rB_hat_2 = [0 -rB[3,2] rB[2,2]; rB[3,2] 0 -rB[1,2]; -rB[2,2] rB[1,2] 0];
    dθ = c .* n
    dθ = dθ * rs_hat;  

    dβ = zeros(numDiodes, 3); # dy/dβ

    dc = n * rs;     # dy/dc  

    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    dα = c .* (ndα * rs);     # dy/dα 

    ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term)
    dϵ = c .* (ndϵ * rs);    # dy/dϵ 


    Cr = [rB_hat_1 zeros(3,3) zeros(3,3*numDiodes);     # Jacobian for gyro measurement portion, [6, 6+3i]
          rB_hat_2 zeros(3,3) zeros(3,3*numDiodes)];

    Ci = [dθ dβ Diagonal(dc) Diagonal(dα) Diagonal(dϵ)]  # Jacobian for current measurement portion, [i, 6+3i]
    C = [Cr; Ci];

    I_meas = c .* (n * rs) .+ 0; # Measured current, no albedo      

    y = [rB[:]; I_meas[:]]; 

    return y
end

function measurement(x, rN, numDiodes, eclipse)
    # Generates the "measured" body vectors for a set of 
    #      newtonian (inertial) vectors using quaternion q
    #      (What our measurement would be given our estimated attitude)
    # Inputs:
    #   - q: quaternion representing current orientation.
    #           Scalar first (10, \vec{q})                             |  [4 x 1]
    #   - rN: Pair of unit vectors in the newtonian (inertial) frame   |  [3 x 2]
    #           (First vector is assumed to be sun vector)
    # Outputs:
    #   - y: Pair of unit vectors in body frame corresponding to rN & q
    #           (Flattened)                                            |  [6 x 1]
    #   - C: Jacobian of y with respect to x
    #           dy/dx = [dy/dPhi; dy/dBeta] = [dy/dPhi; 0]             |  [6 x 6]

    x = x[:]
    q = x[1:4]
    β = x[5:7]


    # DCM from quaternion (flipped)
    B_Q_N = dcm_from_q(q)';                           

    # this is what the measurement would be given our estimated attitude
    rB = B_Q_N*rN;

    rs = B_Q_N * rN[:,1]; # Assumes that rN[:,1] is the sun vector ##################
    rs_hat = [0     -rs[3]   rs[2];   rs[3]    0    -rs[1];  -rs[2]  rs[1]     0];

    c = x[8:(7+numDiodes)]  
    α = x[(8+numDiodes):(7+2*numDiodes)]
    ϵ = x[(8+2*numDiodes):end]

    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];

    rB_hat_1 = [0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    rB_hat_2 = [0 -rB[3,2] rB[2,2]; rB[3,2] 0 -rB[1,2]; -rB[2,2] rB[1,2] 0];
    dθ = c .* n
    dθ = dθ * rs_hat;  

    dβ = zeros(numDiodes, 3); # dy/dβ

    dc = n * rs;     # dy/dc  

    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    dα = c .* (ndα * rs);    # dy/dα 

    ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term)
    dϵ = c .* (ndϵ * rs);    # dy/dϵ 


    Cg = [rB_hat_1 zeros(3,3) zeros(3,3*numDiodes);      # Jacobian for gyro measurement portion, [6, 6+3i]
          rB_hat_2 zeros(3,3) zeros(3,3*numDiodes)];

    Ci = [dθ dβ Diagonal(dc) Diagonal(dα) Diagonal(dϵ)]  # Jacobian for current measurement portion, [i, 6+3i]


    I_meas = c .* (n * rs) .+ 0; # Measured current, no albedo    
    

    I_meas[I_meas .≤ 0] .= 0  # Photodiodes don't generate negative current
    Ci[I_meas .≤ 0, :] .= 0    #  ^ Adjust the jacobian accordingly (note that technically it would change for slight angle changes near the boundary, but we ignore that)

    I_meas = I_meas .* eclipse;
    Ci = Ci .* eclipse;

    C = [Cg; Ci];
    y = [rB[:]; I_meas[:]]; 

    return y, C
end