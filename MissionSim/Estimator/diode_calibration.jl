####################################################################
#               DIODE CALIBRATION + MEKF                           #
####################################################################
# Make sketchy globals into a struct (same with mag calib)

using Infiltrator

mutable struct DIODE_CALIB 
    albedo::ALBEDO        
    sat_state     # x (Full)
    # covariance    # P
    inertial_vecs # rᴵ
    ang_vel       # ω
    body_vecs     # rᴮ
    current_meas  # y 
    W   # noise 1, 2
    V
    dt 
    time 
    num_diodes 
    pos
    first_pass # bool flag
end

function estimate_vals(sat::SATELLITE, data::DIODE_CALIB)
    """
        (Summary)

        Arguments:
            - sat ESTIMATE

        Returns:
            -
    """
    i = data.num_diodes
    data.time = (data.dt + data.time)  # - dt? Need offset of one?
    new_state, new_covariance = mekf(data.sat_state, sat.covariance, data.W, data.V,   # CHANGED data. -> sat.
                                        data.inertial_vecs, data.ang_vel, data.body_vecs, 
                                        data.current_meas, data.num_diodes, data.pos, data.dt, 
                                        data.time, data.albedo)

    i = data.num_diodes
    scale_factors = new_state[8:(7+i)]
    azi_angles    = new_state[(8+i):(7+2*i)]
    elev_angles   = new_state[(8+2*i):(end)]

    diodes = DIODES(scale_factors, azi_angles, elev_angles)
    sat.diodes = diodes

    data.sat_state = new_state 
    # data.covariance = new_covariance 

    sat.state = new_state[1:7] # Only tracking the non-calibration states here
    sat.covariance = new_covariance

    return sat, data 
end


function initialize(albedo, state, system) 
        println("Need to not use state in diode init! (or sat_state in TRIAD)")
        sat_state = [state[7:10][:]; zeros(21)[:]]
        d = DIODE_CALIB(albedo,
                        sat_state, #sat_state, # NOT the same as x0
                        # 0.0, # Just empty stuff to be filled in 
                        0.0, # UPDATE with rᴵ
                        state[11:13], # angular velocity 
                        0.0, # UPDATE with rᴮ
                        0.0, # UPDATE with current_meas 
                        0.0, # Update with W 
                        0.0, # Update with V 
                        system._dt, 
                        system._epc,
                        system._num_diodes,
                        state[1:3], # Position
                        true) # First pass
        return d
end

println("Need to fix covariance in new_diode_calib")
function new_diode_calib(albedo, sens::SENSORS, system, q, sat) 
    """
        (Summary)

        Arguments:
            - 

        Returns:
            -
    """
    sat_state = zeros(7 + 3 * system._num_diodes)
    data = DIODE_CALIB(albedo,
                        sat_state, # NOT the same as x0
                        # 0.0, # UPDATE with covariance  <---- does this need to be reset or passed on through eclipse?
                        0.0, # UPDATE with rᴵ
                        sens.gyro, # angular velocity 
                        0.0, # UPDATE with rᴮ
                        0.0, # UPDATE with current_meas 
                        0.0, # Update with W 
                        0.0, # Update with V 
                        system._dt, 
                        system._epc,
                        system._num_diodes,
                        sens.gps, # Position
                        true) # First pass

    q₀ = q
    β₀ = [0; 0; 0]
    x₀ = [q₀; β₀; sat.diodes.calib_values; sat.diodes.azi_angles; sat.diodes.elev_angles]

    σ_q = (10*pi/180) 
    σ_β = (10*pi/180)
    σ_c = 0.2 #0.15
    σ_α = 1.0 #deg2rad(2.0) #1.0;
    σ_ϵ = 0.3 #deg2rad(2.0) #0.3

    if isnan(sat.covariance[1,1]) # First time
        p = [σ_q * ones(3); σ_β * ones(3); σ_c * ones(data.num_diodes); σ_α*ones(data.num_diodes); σ_ϵ*ones(data.num_diodes)].^2
        P₀ = diagm(p)
    else 
        p = [σ_q * ones(3); σ_β * ones(3)].^2   # Reset Attitude and Bias covariance
        P_new = diagm(p)
        P₀ = sat.covariance 
        P₀[1:6, 1:6] = P_new # Keep calibration covariances
    end

    estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
                        gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
                        velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
                        accel_bias_instability = 6)      # in microG

    Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  # Units are now rad^2/seconds^3...? => (rad/sec)^
    σ_orient = sqrt(Q_gyro);

    Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)   # This is super small
    σ_bias = sqrt(Q_bias)

    Q_diode = 1e-6   # Diode Noise  ##############################3
    σ_cal = Q_diode; σ_azi = Q_diode; σ_ele = Q_diode;
    σ_sunVec = deg2rad(5.0); σ_magVec = deg2rad(5.0); σ_curr = 0.008; #3, 3, 0.005

    W = Diagonal([σ_orient * ones(3); σ_bias * ones(3); σ_cal * ones(data.num_diodes); σ_azi * ones(data.num_diodes); σ_ele * ones(data.num_diodes)])
    V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)])
    
    data.sat_state = x₀
    sat.covariance = P₀
    data.W = W 
    data.V = V

    data.first_pass = false 

    return data
end

    function mekf(x, P, W, V, rᴵ, w, rᴮ,  y, _num_diodes, pos, dt, time, alb::ALBEDO)
        """
            (Summary)

            Arguments:
                - 

            Returns:
                -
        """
        # Runs a single step of a multiplicative extended Kalman filter. 
        #   (Note that the code currently assumes no correlation at all in noise matrix)
        #   returns the next x and P values 
        if sum(abs.(rᴵ[1,:])) < 0.01
            println("Eclipse!") # Shouldn't ever happen
            eclipse = true  
        else
            eclipse = false
        end

        # Predict x, P
        # NORMALIZE ALL VECTORS TO MAKE THEM UNIT 

        sᴵ = rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
        Bᴵ = rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)
        𝐬ᴮ = rᴮ[1,:];  𝐬ᴮ = 𝐬ᴮ / norm(𝐬ᴮ)                    #### CURRENTLY ONLY HAVE UNIT, doing this to be redundantly safe
        Bᴮ = rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

        x_p, A = prediction(x, w, dt, _num_diodes); # State prediction
        P_p = A*P*A' + W; # Covariance prediction 

        # Measurement
        z = [] 
        C = Array{Float64}(undef, 0, size(P_p,1))
        Vs = Array{Float64}(undef)  ## ASSUMES NO CORRELATION AT ALL !!


        # of statements no longer needed!
        if true            #  Magnetic Field measurement
            V_mag = V[1:3, 1:3];
            yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes)   
            z_mag = 𝐁ᴮ - yp_mag 
            z = z_mag[:]
            C = C_mag
            Vs = [diag(V_mag)[:]...]
        end
        
        if !eclipse    #  Diode Current Measurement
            V_cur = V[4:end, 4:end]
            yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ, _num_diodes, 1.0, pos, time, alb) # If !eclipse, set ν = 1.0
            z_cur = y - yp_cur 
            z = [z[:]; z_cur[:]]
            C = [C; C_cur]
            Vs = [Vs; diag(V_cur)[:]]  
        end

        # Innovation   
        # Vs = Vs[2:end] # Get rid of the initial 0 term  
        Vk = Diagonal(Vs) 
        S = C*P_p*C' + Vk; 

        if rank(S) != size(S, 1)
            @infiltrate 
        end


        # Kalman Gain
        L = P_p * C' * S^(-1); 

        # Update
        dx = L*z;          
        dPhi = dx[1:3]; 
        drest = dx[4:end]

        theta_temp = (norm(dPhi));
        rTemp = dPhi / theta_temp; 
        
        dq = [rTemp*sin(theta_temp/2); cos(theta_temp/2)];

        x_next = deepcopy(x)
        x_next[1:4] = qmult(x_p[1:4], dq); 
        x_next[5:end] = x_p[5:end] + drest;
        
        P_next = (I(size(P,1)) - L*C) * P_p * (I(size(P,1)) - L*C)' + L*Vk*L';  

        return x_next, P_next
    end

        function mag_measurement(x, 𝐁ᴵ, i)
            """
                (Summary)

                Arguments:
                    - 

                Returns:
                    -
            """
            # Generates the "measured" body vector from the magnetic field vector
            #      (What our measurement would be given our estimated attitude)
            # Inputs:
            #   - x: Current state of the satellite 
            #           ([q⃗, q0] β⃗ C⃗ α⃗ ϵ⃗)                                            |  [6 + 3i,]
            #   - 𝐁ᴵ: Unit magnetic field vector in the newtonian (inertial) frame   |  [3,]
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

            ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (flipped)    
            𝐁ᴮ = ᴮQᴵ*𝐁ᴵ;     # this is what the measurement would be given our estimated attitude

            B̂ᴮ = hat(𝐁ᴮ); # Hat as in skew-symmetric, not unit or estimate 

            ∂θ = B̂ᴮ 
            ∂β = zeros(3, 3)
            ∂C = zeros(3, i)
            ∂α = zeros(3, i)
            ∂ϵ = zeros(3, i)    

            H = [∂θ ∂β ∂C ∂α ∂ϵ]; # [3 x 6 + 3i]
            y = 𝐁ᴮ[:]             # [3 x 1]

            return y, H
        end

        function current_measurement(x, 𝐬ᴵ, i, ecl, pos, time, alb::ALBEDO)
            """
                (Summary)

                Arguments:
                    - 

                Returns:
                    -
            """
            # Generates the "measured" current values from the sun vector
            #      (What our measurement would be given our sun vector)
            # Inputs:
            #   - x: Current state of the satellite 
            #           ([q⃗, q0] β⃗ C⃗ α⃗ ϵ⃗)                                    |  [6 + 3i,]
            #   - 𝐬ᴵ: Unit sun vector in the newtonian (inertial) frame      |  [3,]
            #   - i: Number of diodes                                        |  (scalar)
            #   - ecl: Eclipse factor η ∈ [0, 1]                             |  (scalar)
            # Outputs:
            #   - y: Current measurements corresponding to 𝐬ᴵ & q            |  [i]
            #   - H: Jacobian of y with respect to x
            #           dy/dx = [dy/dPhi; dy/dBeta; ...]                     |  [i x 6 + 3i]

            x = x[:]
            q = x[1:4]
            β = x[5:7]
            c = x[8:(7+i)]  
            α = x[(8+i):(7+2*i)]
            ϵ = x[(8+2*i):end]

            ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (transposed to get I -> B)
            sᴮ = ᴮQᴵ * 𝐬ᴵ
        
            ŝᴮ= hat(sᴮ);  # use ̌  \check?
            n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
            
            ∂θ = (c .* n) * ŝᴮ; # [i x 3]
            ∂β = zeros(i, 3);       # [i x 3]  
            ∂C = n * sᴮ;            # [i x 1]
        
            ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
            ∂α = c .* (ndα * sᴮ);   # [i x 1]
        
            ndϵ = [(-sin.(ϵ).*cos.(α)) (-1*(sin.(ϵ).*sin.(α))) cos.(ϵ)]; # (With negative middle term)
            ∂ϵ = c .* (ndϵ * sᴮ);   # [i x 1]  
        
            H = [∂θ ∂β Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]
        
        
            I_meas = c .* (n * sᴮ) .+ 0; # Measured current, ALBEDO added in later
        
        
            # ADD IN ALBEDO
            sᴵ_unscaled = sun_position(time) - pos;
            ecl = eclipse_conical(-pos, sᴵ_unscaled) ####### NEED TO FIX TO +pos when updated
            if ecl > 0.98
                ecl = 1.0
            else
                ecl = 0.0
            end
        
            albedo_matrix, ignore = albedo(pos, sᴵ_unscaled, alb.refl)
        
            for j = 1:i
                surface_normal = [cos(ϵ[j])*cos(α[j]) cos(ϵ[j])*sin(α[j]) sin(ϵ[j])]     # Photodiode surface normal 
        
                diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, pos)
        
                diode_albedo = c[j] * diode_albedo / _E_am0;
        
                I_meas[j] = I_meas[j] + diode_albedo
            end
            
        
            # Account for eclipses
            I_meas *= ecl
            I_meas[I_meas .< 0] .= 0  # Photodiodes don't generate negative current
            H[I_meas .≤ 0, :] .= 0    # ^ To match the above
            y = I_meas[:]             # [i x 1]        
            return y, H
        end

        function prediction(xk, w, dt, numDiodes)
            """
                (Summary)

                Arguments:
                    - 

                Returns:
                    -
            """
            # Predicts next state using current state and angular velocity
            #
            # Arguments:
            #   - xk: Current state [q β c α ϵ]                  | [7 + 3i,]
            #   - w: Current angular velocity                    | [3,]
            #   - dt: Time step                                  | Scalar
            #   - numDiodes: Number of photodiodes being used    | Scalar
            #
            # Returns:
            #   - xn: Predicted next state [q β c α ϵ]                     | [7 + 3i,]
            #   - A: State Jacobian with respect to state    
            #           (note that quaternions are replaced with 3 param)  | [(6 + 3i) x (6 + 3i)]
        
            q = xk[1:4]; # Quaternion portion
            b = xk[5:7]; # Bias portion
        
            γ = w-b;     # Adjusted angular velocity (w - biases)
            nγ = norm(γ)
        
            theta = (nγ*dt);  
            r = γ/nγ;  # Make unit
        
            qp = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
            
            skew = -hat(γ)

            R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues (for matrix exponential?)
        
            A = [R -dt*I(3); zeros(3,3) I(3)]; # Jacobian of f(x)
        
            c = xk[8:(7+numDiodes)]
            α = xk[(8+numDiodes):(7 + 2*numDiodes)]
            ϵ = xk[(8 + 2*numDiodes):end]
        
            A = [A                           zeros(6, 3 * numDiodes); 
                zeros(3 * numDiodes, 6)     I(3*numDiodes)];
        
            xn = [qp; b; c; α; ϵ]; # x at next step
        
            return xn, A
        end

     