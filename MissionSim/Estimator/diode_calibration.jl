####################################################################
#               DIODE CALIBRATION + MEKF                           #
####################################################################
# Make sketchy globals into a struct (same with mag calib)

using Infiltrator


mutable struct DIODE_CALIB 
    albedo::ALBEDO        
    sat_state     # x
    covariance    # P
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

    # STATE here != state of simulator !!!!
function estimate_vals(sat::SATELLITE, data::DIODE_CALIB)
    # data:
    #   state, covariance, W, V, Newtonian vectors, w, body vectors, currents, 
    #           (dt, epc, num_diodes?)
    # state: [q β c α ϵ]

    𝐬ᴵ = data.inertial_vecs[1, :] / norm(data.inertial_vecs[1, :])      
    𝐁ᴵ = data.inertial_vecs[2, :] / norm(data.inertial_vecs[2, :])      
    𝐬ᴮ = data.body_vecs[1, :] / norm(data.body_vecs[1, :])              # Should already be unit 
    Bᴮ = data.body_vecs[2, :] / norm(data.body_vecs[2, :])              

    if data.first_pass
        q̂, R₀ = triad(𝐬ᴵ, 𝐁ᴵ, 𝐬ᴮ, Bᴮ)

        # q = data.sat_state[1:4]
        # println("Using state knowledge to init q₀")
        q₀ = q̂


        β₀ = [0; 0; 0]
        x₀ = [q₀; β₀; sat.diodes.calib_values; sat.diodes.azi_angles; sat.diodes.elev_angles]

        σ_q = (10*pi/180)
        σ_β = (10*pi/180)
        σ_c = 0.2; σ_α = 1.0; σ_ϵ = 0.3
        p = [σ_q * ones(3); σ_β * ones(3); σ_c * ones(data.num_diodes); σ_α*ones(data.num_diodes); σ_ϵ*ones(data.num_diodes)].^2
        P₀ = diagm(p)

        # TRY !!
        # P₀ = (10 * pi / 180)^2 * I(6 + 3*6)

        estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
                            gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
                            velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
                            accel_bias_instability = 6)      # in microG

        Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  # Units are now rad^2/seconds^3...?
        σ_orient = sqrt(Q_gyro);

        Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)
        σ_bias = sqrt(Q_bias)

        Q_diode = 1e-6 # Diode Noise 

        σ_cal = Q_diode; σ_azi = Q_diode; σ_ele = Q_diode;
        σ_sunVec = deg2rad(5.0); σ_magVec = deg2rad(5.0); σ_curr = 0.008;

        W = Diagonal([σ_orient * ones(3); σ_bias * ones(3); σ_cal * ones(data.num_diodes); σ_azi * ones(data.num_diodes); σ_ele * ones(data.num_diodes)])
        # V = Diagonal([σ_sunVec * ones(3); σ_magVec * ones(3); σ_curr * ones(data.num_diodes)])
        V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)])
        
        data.sat_state = x₀
        data.covariance = P₀
        data.W = W 
        data.V = V

        data.first_pass = false 
    else
        # println("In MEKF")
        data.time = data.dt + data.time
    end
    
        new_state, new_covariance = mekf(data.sat_state, data.covariance, data.W, data.V,
                                            data.inertial_vecs, data.ang_vel, data.body_vecs, 
                                            data.current_meas, data.num_diodes, data.pos, data.dt, 
                                            data.time, data.albedo)

        i = data.num_diodes
        scale_factors = wrap(new_state[8:(7+i)])
        azi_angles    = wrap(new_state[(8+i):(7+2*i)])
        elev_angles   = wrap(new_state[(8+2*i):(end)])

        diodes = DIODES(scale_factors, azi_angles, elev_angles)
        sat.diodes = diodes

        data.sat_state = new_state 
        data.covariance = new_covariance 
    # end

    return sat, data, false
end

    function wrap(angles)
        # angles[angles .> (2 * pi)] .= angles[angles .> (2 * pi)] .- (2 * pi)
        # angles[angles .< 0] .= angles[angles .< 0] .+ (2 * pi)
        return angles
    end

function initialize(albedo, state, system) 
    println("Need to not use state in diode init! (or sat_state in TRIAD)")
    sat_state = [state[7:10][:]; zeros(21)[:]]
    d = DIODE_CALIB(albedo,
                    sat_state, #sat_state, # NOT the same as x0
                    0.0, # Just empty stuff to be filled in 
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

function new_diode_calib(albedo, sens::SENSORS, system) 
    d = DIODE_CALIB(albedo,
                    0.0, # NOT the same as x0
                    0.0, # Just empty stuff to be filled in 
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
    return d
end

    function mekf(x, P, W, V, rᴵ, w, rᴮ,  y, _num_diodes, pos, dt, time, alb::ALBEDO)
        # Runs a single step of a multiplicative extended Kalman filter. 
        #   (Note that the code currently assumes no correlation at all in noise matrix)
        #   returns the next x and P values 

        if sum(abs.(rᴵ[1,:])) < 0.01
            eclipse = true  
        else
            eclipse = false
        end

        # Predict x, P
        # NORMALIZE ALL VECTORS TO MAKE THEM UNIT 
        # rᴵ not unit, rᴮ[1,:] already unit

        sᴵ = rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
        Bᴵ = rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)
        𝐬ᴮ = rᴮ[1,:]                      #### CURRENTLY ONLY HAVE UNIT
        Bᴮ = rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

        x_p, A = prediction(x, w, dt, _num_diodes); # State prediction
        P_p = A*P*A' + W; # Covariance prediction 

        # Measurement
        z = [] 
        C = Array{Float64}(undef, 0, size(P_p,1))
        Vs = Array{Float64}(undef)  ## ASSUMES NO CORRELATION AT ALL !!

        # if false #!eclipse  #  Sun Vector measurement
            #     V_sun = V[1:3, 1:3]; 
            #     yp_sun, C_sun = sun_measurement(x_p, rᴵ[1, :], _num_diodes) 
            #     z_sun = rᴮ[1,:]  - yp_sun
            #     z = [z[:]; z_sun[:]]
            #     C = [C; C_sun]
            #     Vs = [Vs[:]; diag(V_sun)[:]]  # Does this need Vs[:]?
        # end

        if true            #  Magnetic Field measurement
            V_mag = V[1:3, 1:3];
            yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes)   
            z_mag = 𝐁ᴮ - yp_mag 
            z = z_mag[:]
            C = C_mag
            Vs = [diag(V_mag)[:]]
        end
        

        if !eclipse    #  Diode Current Measurement
            V_cur = V[4:end, 4:end]
            yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ, _num_diodes, 1.0, pos, time, alb) # If !eclipse, set ν = 1.0
            z_cur = y - yp_cur 
            z = [z[:]; z_cur[:]]
            C = [C; C_cur]
            Vs = [Vs...; diag(V_cur)[:]]  
        end


        # Innovation   
        # Vs = Vs[2:end] # Get rid of the initial 0 term  
        Vk = Diagonal(Vs) 
        S = C*P_p*C' + Vk;  

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

            B_Q_N = dcm_from_q(q)'; # DCM from quaternion (flipped)   TODO why flipped?    
            sB = B_Q_N*sN;  # this is what the sun measurement would be given our estimated attitude

            sB_hat = hat(sB); 
            ∂θ = sB_hat 
            ∂β = zeros(3, 3)
            ∂C = zeros(3, i)
            ∂α = zeros(3, i)
            ∂ϵ = zeros(3, i)

            H = [∂θ ∂β ∂C ∂α ∂ϵ]; # [3 x 6 + 3i]
            y = sB[:]             # [3 x 1]

            return y, H
        end

        function mag_measurement(x, 𝐁ᴵ, i)
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

            ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (transposed)
            sᴮ = ᴮQᴵ * 𝐬ᴵ
        
            ŝᴮ= hat(sᴮ); 
            n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
            
            ∂θ = (c .* n) * ŝᴮ; # [i x 3]
            ∂β = zeros(i, 3);       # [i x 3]
            ∂C = n * sᴮ;            # [i x 1]
        
            ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
            ∂α = c .* (ndα * sᴮ);   # [i x 1]
        
            ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term)
            ∂ϵ = c .* (ndϵ * sᴮ);   # [i x 1]  
        
            H = [∂θ ∂β Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]
        
        
            I_meas = c .* (n * sᴮ) .+ 0; # Measured current, ALBEDO added in later
        
        
            # ADD IN ALBEDO
            sᴵ_unscaled = sun_position(time) - pos;
        
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

        function triad(rN1,rN2,rB1,rB2)
            # Method for estimating the rotation matrix between two reference frames 
            #   Relies only on a single pair of vectors in each frame
            # Inputs: 
            #   - rN1, rN2: Pair of vectors in the Newtonian (inertial) frame     | [3,]
            #   - rB1, rB2: Corresponding pair of vectors in body frame           | [3,]    # 
            # Outputs:
            #   - R: A directed cosine matrix (DCM) representing the rotation     | [3 x 3]
            #           between the two frames 
            #   - q: A quaternion (scalar last) representing the rotation         | [4,]
            #           between the two frames
        
            𝐫₁ᴵ = rN1 / norm(rN1)
            𝐫₂ᴵ = rN2 / norm(rN2)
            𝐫₁ᴮ = rB1 / norm(rB1)
            𝐫₂ᴮ = rB2 / norm(rB2)

            t₁ᴵ = 𝐫₁ᴵ
            t₂ᴵ = cross(𝐫₁ᴵ, 𝐫₂ᴵ)/norm(cross(𝐫₁ᴵ,𝐫₂ᴵ));
            t₃ᴵ = cross(t₁ᴵ, t₂ᴵ)/norm(cross(t₁ᴵ,t₂ᴵ));

            Tᴵ = [t₁ᴵ[:] t₂ᴵ[:] t₃ᴵ[:]]

            t₁ᴮ = 𝐫₁ᴮ
            t₂ᴮ = cross(𝐫₁ᴮ, 𝐫₂ᴮ)/norm(cross(𝐫₁ᴮ,𝐫₂ᴮ));
            t₃ᴮ = cross(t₁ᴮ, t₂ᴮ)/norm(cross(t₁ᴮ,t₂ᴮ));

            Tᴮ = [t₁ᴮ[:] t₂ᴮ[:] t₃ᴮ[:]]

            R = Tᴵ * (Tᴮ')

            # q = q_from_DCM(R);
        
            # return q, R

        
            rN1 = rN1 / norm(rN1)
            rN2 = rN2 / norm(rN2)
            rB1 /= norm(rB1)
            rB2 /= norm(rB2)
            # This method avoids taking an inverse by generating an orthogonal matrix
            tN1 = rN1;
            tN2 = cross(rN1,rN2)/norm(cross(rN1,rN2));
            tN3 = cross(tN1,tN2)/norm(cross(tN1,tN2));  
        
            nT = [tN1[:] tN2[:] tN3[:]];
        
        
            tB1 = rB1;
            tB2 = cross(rB1,rB2)/norm(cross(rB1,rB2));
            tB3 = cross(tB1,tB2)/norm(cross(tB1,tB2))
        
            bT = [tB1[:] tB2[:] tB3[:]];
        
            # DCM
            R = nT*(bT');

        
            # QUATERNION
            q = q_from_DCM(R);
        
            return q, R
        end

    