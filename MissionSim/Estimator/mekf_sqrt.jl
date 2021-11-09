####################################################################
#                      STANDARD MEKF                               #
####################################################################


mutable struct MEKF_DATA 
    albedo::ALBEDO  #  Albedo Struct
    sat_state       #  Current state + calibration estimates [(q⃗, q₀) β C α ϵ]  |  [7 + 3i,]
    # covariance 
    inertial_vecs   #  Reference vectors in inertial frame 
    body_vecs       #  Reference vectors in body frame
    ang_vel         #  Angular velocity      
    current_meas    #  Photodiode current measurements       
    pos             #  Current satellite position estimate
    W               #  Process noise matrix   
    V               #  Measurement noise matrix s
    dt              #  Simulation time step           
    time            #  Current time (as an epoch) 
    num_diodes      #  Number of photodiodes (i)  
    mag_calib_matrix # T 
end

# DOES NOT update calibration parameters
function estimate_vals(sat::SATELLITE, data::MEKF_DATA)
    """
        Calls a single iteration of a multiplicative extended Kalman filter 
                and updates satellite attitude and bias estimates, BUT DOES NOT update calibration parameters

        Arguments:
            - sat:  Struct containing current estimate of satellite values      |  SATELLITE 
            - data: Struct containing necessary information to run MEKF         |  -----

        Returns:
            - sat:  Updated satellite struct containing new estimates           |  SATELLITE
            - data: Updated struct containing necessary MEKF information        |  -----
    """  
    data.time += data.dt
    x = data.sat_state[1:10]
    c, α, ϵ = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles
    
    new_state, new_covariance = mekf_sqrt(x, c, α, ϵ, sat.covariance[1:6, 1:6], data.W, data.V,
                                        data.inertial_vecs, data.body_vecs, data.ang_vel, 
                                        data.current_meas, data.num_diodes, data.pos, data.dt, 
                                        data.time, data.albedo, data.mag_calib_matrix)

    data.sat_state = new_state 
    # data.covariance = new_covariance 

    sat.state = new_state[1:10] # Only tracking the non-calibration states
    sat.covariance[1:9, 1:9] = new_covariance
    return sat, data 
end

function new_mekf_data(data::DIODE_CALIB)
    """
        Generates a new MEKF_DATA struct containing necessary information to estimate attitude and gyro bias.
            Uses an existing DIODE_CALIB struct, meaning this is called when the diodes are calibrated while in sunlight

        Arguments:
            - data:  DIODE_CALIB struct containing necessary information to initialize the vanilla MEKF    | DIODE_CALIB

        Returns:
            - mekf:  MEKF_DATA struct       | MEKF_DATA
    """    
    
    mekf = MEKF_DATA(data.albedo, data.sat_state[1:7], # data.covariance[1:6, 1:6], 
                        data.inertial_vecs,  data.body_vecs, data.ang_vel,
                        data.current_meas, data.pos, data.W[1:6, 1:6], data.V, 
                        data.dt, data.time, data.num_diodes, data.mag_calib_matrix)

    return mekf
end 

function new_mekf_data(alb::ALBEDO, sens::SENSORS, system, q, sat)
    """
        Generates a new MEKF_DATA struct containing necessary information to estimate attitude and gyro bias.
            This version is called when coming out of an eclipse, and must initialize all values 

        Arguments:
            - alb: ALBEDO struct containing refl and cell centers         |  ALBEDO
            - sens: SENSORS struct containing sensor measurements         |  SENSORS 
            - system: Set of system constants needed for Simulation       |  List
            - q: Current estimate of satellite attitude                   |  [4,] ([q⃗; q₀] -> Scalar last)
            - sat: Struct containing current satellite estimate           |  SATELLITE

        Returns:
            - mekf:  MEKF_DATA struct                                     | MEKF_DATA
    """  
    sat_state = zeros(10)
    data = MEKF_DATA(alb,        # ALBEDO 
                        sat_state,    # Satellite state ([(q⃗, q₀) β C α ϵ], NOT the same as state in simulation)
                        # 0.0,        # UPDATE with covariance  <---- does this need to be reset or passed on through eclipse?
                        0.0,          # Ref vectors in inertial frame (rᴵ)
                        0.0,          # Ref vectors in body frame (rᴮ)
                        sens.gyro,    # Angular velocity 
                        0.0,          # Photodiode current measurements 
                        sens.gps,     # Position
                        0.0,          # Update with W 
                        0.0,          # Update with V 
                        system._dt,   # Time step
                        system._epc,  # Initial time
                        system._num_diodes,  # Number of photodiodes on the satellite 
                        generate_mag_calib_matrix(sat))  

    β₀_gyro = sat.state[5:7] # Continue with last estimate of bias
    β₀_mag  = sat.state[8:10]
    x₀ = [q; β₀_gyro; β₀_mag] 

    # σ_q = (10*pi/180) 
    # σ_β = (10*pi/180)


    if isnan(sat.covariance[1,1]) # Setting up covariance from scratch
        p = [σ_q * ones(3); σ_βgyro * ones(3); σ_βmag * ones(3); zeros(18)].^2  # in theory, if we are starting the MEKF from scratch we know the calibration parameters...
        sat.covariance = (diagm(p))
    end
    p = [σ_q * ones(3); σ_βgyro * ones(3); σ_βmag * ones(3)].^2   # Reset Attitude and Bias covariance
    P₀ = diagm(p)

    # ##### VERIFY THIS SECTION ##########  -> IN config file?
    # estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
    #                     gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
    #                     velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
    #                     accel_bias_instability = 6)      # in microG

    # Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  # Units are now rad^2/seconds^3...? => (rad/sec)^
    # σ_orient = sqrt(Q_gyro);

    # Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)   # This is super small
    # σ_bias = sqrt(Q_bias)

    # σ_sunVec = deg2rad(3.0); σ_magVec = deg2rad(3.0); σ_curr = 0.05; #3, 3, 0.005
    # #######################################

    W = Diagonal([σ_orient * ones(3); σ_bias_gyro * ones(3); σ_bias_mag * ones(3)]).^2
    V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)]).^2
    
    data.sat_state = x₀
    # sat.covariance = P₀   #### Should this be sat.covariance[1:6, 1:6] = P₀?
    sat.covariance[1:9, 1:9] = chol(P₀)
    # sat.covariance = chol(P₀[1:6, 1:6])
    data.W = W[1:9, 1:9] # W
    data.V = V
            
    return data
end

function mag_measurement(x, 𝐁ᴵ, T)
    """
        Generates what we would expect our measured magnetic field vector would be in the body frame 
            given our current estimated attitude

        Arguments:
            - x:  Current state of the satellite [(q⃗, q₀) β C α ϵ]             | [7 + 3i,] 
            - 𝐁ᴵ: Unit mag field vector in the newtonian (inertial) frame      |  [3,]

        Returns:
            - y: Unit mag field vector estimate corresponding to state         |  [3,]
            - H: Jacobian of y with respect to x                               |  [3 x 6 + 3i]
                        dy/dx = [dy/dϕ; dy/dβ; ...]                             
    """

    x = x[:]
    q = x[1:4]
    β_gyro = x[5:7]
    β_mag  = x[8:10]

    ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (flipped)    
    𝐁ᴮ = ᴮQᴵ*𝐁ᴵ;     # this is what the measurement would be given our estimated attitude

    B̌ᴮ = hat(𝐁ᴮ);    # Skew-symmetric matrix

    ∂θ = B̌ᴮ 
    ∂βgyro = zeros(3, 3)
    ∂βmag  = T

    H = [∂θ ∂βgyro ∂βmag]; # [3 x 6]
    y = 𝐁ᴮ[:]    # [3 x 1]

    return y, H
end

function current_measurement(x, c, α, ϵ, 𝐬ᴵ, i, pos, time, alb::ALBEDO)
    """
        Generates what we would expect our measured current values to be given our estimate of state

        Arguments:
            - x:  Current state of the satellite [(q⃗, q₀) β C α ϵ]             | [7 + 3i,] 
            - c: Best estimate of diode calibration values                     |  [i,]
            - α: Best estimate of diode azimuth angles                         |  [i,]
            - ϵ: Best estimate of diode elevation angles                       |  [i,]    
            - 𝐬ᴵ: Unit sun vector in the newtonian (inertial) frame            |  [3,]
            - i:  Number of photodiodes                                        |  Int
            - pos: Current satellite position estimate                         |  [3,]
            - time: Current time (as an epoch)                                 |  Epoch 
            - alb: ALBEDO struct containing refl and cell centers              |  ALBEDO

        Returns:
            - y: Current measurements corresponding to 𝐬ᴵ & q (non-negative)   |  [i,]
            - H: Jacobian of y with respect to x                               |  [i x 6 + 3i]
                        dy/dx = [dy/dϕ; dy/dβ; ...]                             
    """
    x = x[:]
    q = x[1:4]
    β = x[5:7]

    ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (transposed to get I -> B)
    sᴮ = ᴮQᴵ * 𝐬ᴵ

    šᴮ= hat(sᴮ);  # Skew-symmetric form
    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
    
    ∂θ = (c .* n) * šᴮ;  # [i x 3]
    ∂βgyro = zeros(i, 3);    # [i x 3]  
    ∂βmag  = zeros(i, 3)
    H = [∂θ ∂βgyro ∂βmag]          # [i x 6]                                           

    I_meas = c .* (n * sᴮ); # Measured current, ALBEDO added in later

    # ADD IN ALBEDO
    sᴵ_unscaled = sun_position(time) - pos;
    # ecl = eclipse_conical(-pos, sᴵ_unscaled) ####### NEED TO FIX TO +pos when updated
    # ecl = (ecl > 0.98) ? 1.0 : 0.0

    albedo_matrix, _ = albedo(pos, sᴵ_unscaled, alb.refl)

    for j = 1:i
        surface_normal = [cos(ϵ[j])*cos(α[j]) cos(ϵ[j])*sin(α[j]) sin(ϵ[j])]     # Photodiode surface normal 
        diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, pos)
        diode_albedo = c[j] * diode_albedo / _E_am0;
        I_meas[j] = I_meas[j] + diode_albedo
    end
    

    # Account for eclipses
    # I_meas *= ecl
    I_meas[I_meas .< 0] .= 0  # Photodiodes don't generate negative current
    H[I_meas .≤ 0, :] .= 0    # ^ To match the above
    y = I_meas[:]             # [i x 1]        
    return y, H
end

function prediction(x, w, dt)
    """
        Predicts the next state and covariance using current state, angular velocity, and time step
            (essentially assumes a small rotation in attitude and that the other states are constant)

        Arguments:
            - x:  Current state of the satellite [(q⃗, q₀) β C α ϵ]              | [7 + 3i,] 
            - w:  Current angular velocity estimate of the satellite            | [3,]
            - dt: Time step of the Simulation                                         |  Int

        Returns:
            - xn: Predicted next state  [(q⃗, q₀) β C α ϵ]                       | [7 + 3i,] 
            - H:  Jacobian of state x with respect to itself                    |  [i x 6 + 3i]
                    (Note that the quaternions have been replaced with 3 param)
                        dx/dx = [dϕ/dϕ; dϕ/dβ; ...]                             
    """
    q = x[1:4]; # Quaternion portion
    βgyro = x[5:7]; # Bias portion
    βmag  = x[8:10]

    γ = w - βgyro;     # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    θ = (nγ*dt);  
    r = γ/nγ;  # Make unit

    qp = qmult(q, [r*sin(θ/2); cos(θ/2)]); 
    
    skew = -hat(γ)

    R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues (for matrix exponential?)

    # A = [R -dt*I(3); zeros(3,3) I(3)]; 
    A = [R         -dt*I(3)   zeros(3,3);  # Jacobian of f(x)
        zeros(3,3) I(3)      zeros(3,3);
        zeros(3,3) zeros(3,3)      I(3)] 

    xn = [qp; βgyro; βmag]  # x at next step

    return xn, A
end

function mekf_sqrt(x, c, α, ϵ, Pchol, W, V, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO, T) 
    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)                
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Prediction
    x_p, A = prediction(x, w, dt) # State prediction
    Pchol_p = qrᵣ([Pchol * A'; chol(Matrix(W))])  # Cholesky factors


    # Measurement
    yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, T)
    z_mag = 𝐁ᴮ - yp_mag 
    yp_cur, C_cur = current_measurement(x_p, c, α, ϵ, 𝐬ᴵ, _num_diodes, pos, time, alb) 
    z_cur = y - yp_cur 

    C = [C_mag; C_cur]
    z = [z_mag[:]; z_cur[:]]

    # Innovation   
    Vk = Diagonal(V) 
    # S = C*P_p*C' + Vk;  
    
    # Kalman Gain
    Pyy_chol = qrᵣ([Pchol_p * C'; chol(Matrix(V))])  
    # L = P_p * C' * S^(-1); 
    L = ((Pchol_p' * Pchol_p * C') / Pyy_chol) / (Pyy_chol')

    # Update
    dx = L*z;    
    dϕ = dx[1:3]; 
    drest = dx[4:end]

    θ = (norm(dϕ));
    rTemp = dϕ / θ; 

    dq = [rTemp*sin(θ/2); cos(θ/2)];

    x_next = deepcopy(x)
    x_next[1:4] = qmult(x_p[1:4], dq);
    x_next[5:end] = x_p[5:7] + drest;

    # T = (I(size(P,1)) - L*C) 
    # P_next = T * P_p * T' + L*Vk*L';  
    Pchol_next = qrᵣ([Pchol_p*(I-L*C)'; chol(Matrix(V))*L'] )

    return x_next, Pchol_next
end