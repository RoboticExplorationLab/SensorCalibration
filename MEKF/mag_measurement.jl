function mag_measurement(x, bN, i)
    # Generates the "measured" body vector from the magnetic field vector
    #      (What our measurement would be given our estimated attitude)
    # Inputs:
    #   - x: Current state of the satellite 
    #           ([q⃗, q0] β⃗ C⃗ α⃗ ϵ⃗)                                            |  [6 + 3i,]
    #   - bN: Unit magnetic field vector in the newtonian (inertial) frame   |  [3,]
    #   - i: Number of diodes                                                |  (scalar)
    # Outputs:
    #   - y: Unit vector in body frame corresponding to bN & q               |  [3,]
    #   - H: Jacobian of y with respect to x
    #           dy/dx = [dy/dPhi; dy/dBeta; ...]                             |  [3 x 6 + 3i]


    x = x[:]
    q = x[1:4]
    β = x[5:7]
    c = x[8:(7+i)]  
    α = x[(8+i):(7+2*i)]
    ϵ = x[(8+2*i):end]

    B_Q_N = dcm_from_q(q)'; # DCM from quaternion (flipped)    
    bB = B_Q_N*bN;     # this is what the measurement would be given our estimated attitude

    bB_hat = hat(bB); 
    ∂θ = bB_hat 
    ∂β = zeros(3, 3)
    ∂C = zeros(3, i)
    ∂α = zeros(3, i)
    ∂ϵ = zeros(3, i)

    H = [∂θ ∂β ∂C ∂α ∂ϵ]; # [3 x 6 + 3i]
    y = bB[:]             # [3 x 1]

    return y, H
end