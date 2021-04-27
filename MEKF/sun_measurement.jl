function sun_measurement(x, sN, i)
    # Generates the "measured" body vector from the sun vector
    #      (What our measurement would be given our estimated attitude)
    # Inputs:
    #   - x: Current state of the satellite 
    #           ([q⃗, q0] β⃗ C⃗ α⃗ ϵ⃗)                                    |  [6 + 3i,]
    #   - sN: Unit sun vector in the newtonian (inertial) frame      |  [3,]
    #   - i: Number of diodes                                        |  (scalar)
    # Outputs:
    #   - y: Unit vector in body frame corresponding to sN & q       |  [3,]
    #   - H: Jacobian of y with respect to x
    #           dy/dx = [dy/dPhi; dy/dBeta; ...]                     |  [3 x 6 + 3i]


    x = x[:]
    q = x[1:4]
    β = x[5:7]
    c = x[8:(7+i)]  
    α = x[(8+i):(7+2*i)]
    ϵ = x[(8+2*i):end]

    B_Q_N = dcm_from_q(q)'; # DCM from quaternion (flipped)    
    sB = B_Q_N*sN;  # this is what the sun measurement would be given our estimated attitude

    sB_hat = hat(sB); #[0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    ∂θ = sB_hat 
    ∂β = zeros(3, 3)
    ∂C = zeros(3, i)
    ∂α = zeros(3, i)
    ∂ϵ = zeros(3, i)

    H = [∂θ ∂β ∂C ∂α ∂ϵ]; # [3 x 6 + 3i]
    y = sB[:]             # [3,]

    return y, H
end