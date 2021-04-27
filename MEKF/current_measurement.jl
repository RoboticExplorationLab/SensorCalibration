function current_measurement(x, sN, i, ecl)
    # Generates the "measured" current values from the sun vector
    #      (What our measurement would be given our sun vector)
    # Inputs:
    #   - x: Current state of the satellite 
    #           ([q⃗, q0] β⃗ C⃗ α⃗ ϵ⃗)                                    |  [6 + 3i,]
    #   - sN: Unit sun vector in the newtonian (inertial) frame      |  [3,]
    #   - i: Number of diodes                                        |  (scalar)
    #   - ecl: Eclipse factor η ∈ [0, 1]                             |  (scalar)
    # Outputs:
    #   - y: Current measurements corresponding to sN & q            |  [i]
    #   - H: Jacobian of y with respect to x
    #           dy/dx = [dy/dPhi; dy/dBeta; ...]                     |  [i x 6 + 3i]


    x = x[:]
    q = x[1:4]
    β = x[5:7]
    c = x[8:(7+i)]  
    α = x[(8+i):(7+2*i)]
    ϵ = x[(8+2*i):end]

    B_Q_N = dcm_from_q(q)'; # DCM from quaternion (flipped)    
    sB = B_Q_N*sN;  # this is what the sun measurement would be given our sun vector

    sB_hat = hat(sB); #[0 -rB[3,1] rB[2,1]; rB[3,1] 0 -rB[1,1]; -rB[2,1] rB[1,1] 0];
    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
    
    ∂θ = (c .* n) * sB_hat; # [i x 3]
    ∂β = zeros(i, 3);       # [i x 3]
    ∂C = n * sB;            # [i x 1]

    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    ∂α = c .* (ndα * sB);   # [i x 1]

    ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term)
    ∂ϵ = c .* (ndϵ * sB);   # [i x 1]  

    H = [∂θ ∂β Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]

    I_meas = c .* (n * sB) .+ 0; # Measured current, NO ALBEDO

    # Account for eclipses
    I_meas *= ecl
    I_meas[I_meas .≤ 0] .= 0  # Photodiodes don't generate negative current
    H[I_meas .≤ 0, :] .= 0    # ^ To match the above
    y = I_meas[:]     # [i,]

    return y, H
end