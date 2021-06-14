function current_measurement(x, sN, i, ecl, pos, time)
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

    sB_hat = hat(sB); 
    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
    
    ∂θ = (c .* n) * sB_hat; # [i x 3]
    ∂β = zeros(i, 3);       # [i x 3]
    ∂C = n * sB;            # [i x 1]

    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    ∂α = c .* (ndα * sB);   # [i x 1]

    ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term)
    ∂ϵ = c .* (ndϵ * sB);   # [i x 1]  

    H = [∂θ ∂β Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]


    I_meas = c .* (n * sB) .+ 0; # Measured current, ALBEDO added in later


    
    #### ACCOUNT FOR ALBEDO... SOMEHOW - do I use pos/epc? Or just use unscaled sN?
    # pos = pos .+ randn(3) * pos * 0.001  # Make it an estimate...?

    sN_unscaled = sun_position(time) - pos;

    albedo_matrix, ignore = albedo(pos, sN_unscaled, refl)

    diode_albedos = zeros(num_diodes)

    for i = 1:num_diodes
        # Below is just the rows of n

        surface_normal = [cos(ϵ[i])*cos(α[i]) cos(ϵ[i])*sin(α[i]) sin(ϵ[i])]     # Photodiode surface normal 

        diode_albedo = get_diode_albedo_local(albedo_matrix, surface_normal, pos)

        diode_albedo = c[i] * diode_albedo / _E_am0;

        I_meas[i] = I_meas[i] + diode_albedo
        diode_albedos[i] = diode_albedo
    end


    #####################################
    

    

    # Account for eclipses
    I_meas *= ecl
    I_meas[I_meas .< 0] .= 0  # Photodiodes don't generate negative current
    H[I_meas .≤ 0, :] .= 0    # ^ To match the above
    y = I_meas[:]     # [i x 1]

    return y, H, diode_albedos[:]
end