####################################################################
#                      STANDARD MEKF                               #
####################################################################

@info "Not sure I am splitting up covariance into the three matrices correctly! (Shouldn't be diagonal...?)"
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
    x = data.sat_state[1:7]
    c, α, ϵ = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles
    
    new_state, new_covariance = mekf(x, c, α, ϵ, sat.covariance, data.W, data.V, 
                                        data.inertial_vecs, data.body_vecs, data.ang_vel, 
                                        data.current_meas, data.num_diodes, data.pos, data.dt, 
                                        data.time, data.albedo)

    data.sat_state = new_state 

    sat.state = new_state[1:7] # Only tracking the non-calibration states
    sat.covariance = new_covariance
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
                        data.dt, data.time, data.num_diodes)

    # mekf = MEKF(data.albedo, data.sat_state[1:7], #.covariance, 
    #             data.inertial_vecs, data.ang_vel, data.body_vecs,
    #             data.current_meas, data.W, data.V, data.dt,
    #             data.time, data.num_diodes, data.pos)
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
    sat_state = zeros(7)
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
                        system._num_diodes)   # Number of photodiodes on the satellite 

    β₀ = sat.state[5:7] # Continue with last estimate of bias
    x₀ = [q; β₀] 

    σ_q = (10*pi/180) 
    σ_β = (10*pi/180)

    if isnan(sat.covariance[1,1]) # Setting up covariance from scratch
        p = [σ_q * ones(3); σ_β * ones(3); zeros(18)].^2
        P₀ = diagm(p)
        sat.covariance = P₀
    else    # Setting up covariance after an eclipse
        p = [σ_q * ones(3); σ_β * ones(3)].^2   # Reset Attitude and Bias covariance
        P₀ = diagm(p)
        sat.covariance[1:6, 1:6] = P₀
    end

    ##### VERIFY THIS SECTION ##########
    estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
                        gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
                        velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
                        accel_bias_instability = 6)      # in microG

    Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  # Units are now rad^2/seconds^3...? => (rad/sec)^
    σ_orient = sqrt(Q_gyro);

    Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)   # This is super small
    σ_bias = sqrt(Q_bias)

    σ_sunVec = deg2rad(5.0); σ_magVec = deg2rad(5.0); σ_curr = 0.008; #3, 3, 0.005
    #######################################

    W = Diagonal([σ_orient * ones(3); σ_bias * ones(3)])

    V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)])
    
    data.sat_state = x₀
    sat.covariance[1:6, 1:6] = P₀   #### Should this be sat.covariance[1:6, 1:6] = P₀?
    data.W = W[1:6, 1:6] 
    data.V = V
            
    return data
end

function mekf(x, c, α, ϵ, Pₐ, W, V, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO)
    """
        Runs a single step of a multiplicative extended Kalman filter 
            (DOES NOT update the calibration parameters)

        Arguments:
            - x: Current state + calibration estimates [(q⃗, q₀) β C α ϵ]  |  [7 + 3i,]
            - c: Best estimate of diode calibration values                |  [i,]
            - α: Best estimate of diode azimuth angles                    |  [i,]
            - ϵ: Best estimate of diode elevation angles                  |  [i,]
            - P: Current covariance matrix for the state                  |  [6 + 3i  x  6 + 3i]
            - W: Process noise matrix                                     |  [6 + 3i  x  6 + 3i]
            - V: Measurement noise matrix                                 |  [3 + i,]
            - rᴵ: Reference vectors in inertial frame                     |  [2 x 3]
            - rᴮ: Reference vectors in body frame                         |  [2 x 3]
            - w: Angular velocity                                         |  [3,]
            - y: Current measurements                                     |  [i,]
            - _num_diodes: Number of photodiodes (i)                      |  Int 
            - pos: Current satellite position estimate                    |  [3,]
            - dt: Simulation time step                                    |  Scalar 
            - time: Current time (as an epoch)                            |  Epoch 
            - alb: ALBEDO struct containing refl and cell centers         |  ALBEDO

        Returns:
            - x_next: Estimate of next state                              |  [7 + 3i,]
            - P_next: Covariance of next state                            |  [6 +3i  x  6 + 3i]
    """      

    P = Pₐ[1:6, 1:6]

    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)
    𝐬ᴮ = @view rᴮ[1,:];  𝐬ᴮ = 𝐬ᴮ / norm(𝐬ᴮ)                    
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Predict x, P
    x_p, Aₐ = prediction(x, w, dt, _num_diodes, false); # State prediction ( WITH EXTRA ARGUMENT TO DISTINGUISH)
    A = Aₐ[1:6, 1:6]

    P_p = A*P*A' + W; # Covariance prediction 

    # Measurement
    yp_mag, C_magₐ = mag_measurement(x_p, c, α, ϵ, 𝐁ᴵ, _num_diodes)
    z_mag = 𝐁ᴮ - yp_mag 
    yp_cur, C_curₐ = current_measurement(x_p, c, α, ϵ, 𝐬ᴵ, _num_diodes, pos, time, alb) 
    z_cur = y - yp_cur 

    Cₐ = [C_magₐ; C_curₐ]
    C = Cₐ[:, 1:6]
    z = [z_mag[:]; z_cur[:]]

    # Innovation   
    Vk = Diagonal(V)
    S = C*P_p*C' + Vk;  
    
    # Kalman Gain
    L = P_p * C' * S^(-1); 

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

    

    Wₐ = [W zeros(6, 18);
            zeros(18, 24) ]

    Kₐ = [L; zeros(18, 9)]

    Pₐ[1:6, 1:6] = P_p

    Iₚ = I(size(Pₐ,1))
    P⁺ = (Iₚ - Kₐ*Cₐ) * Pₐ * (Iₚ - Kₐ*Cₐ) + Kₐ*Vk*Kₐ';
    
    # P_next = (I(size(P,1)) - L*C) * P_p * (I(size(P,1)) - L*C)' + L*Vk*L';  

    return x_next, P⁺
end

function mag_measurement(x, c, α, ϵ, 𝐁ᴵ, i)
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
    β = x[5:7]

    ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (flipped)    
    𝐁ᴮ = ᴮQᴵ*𝐁ᴵ;     # this is what the measurement would be given our estimated attitude

    B̌ᴮ = hat(𝐁ᴮ);    # Skew-symmetric matrix

    ∂θ = B̌ᴮ 
    ∂β = zeros(3, 3)
    ∂C = zeros(3, i)
    ∂α = zeros(3, i)
    ∂ϵ = zeros(3, i)    

    H = [∂θ ∂β ∂C ∂α ∂ϵ]; # [3 x 6 + 3i]
    y = 𝐁ᴮ[:]             # [3 x 1]

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
    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    ndϵ = [(-sin.(ϵ).*cos.(α)) (-1*(sin.(ϵ).*sin.(α))) cos.(ϵ)]; # (With negative middle term, differing from the paper)
    
    ∂θ = (c .* n) * šᴮ;     # [i x 3]
    ∂β = zeros(i, 3);       # [i x 3]  
    ∂C = n * sᴮ;            # [i,]
    ∂α = c .* (ndα * sᴮ);   # [i,]
    ∂ϵ = c .* (ndϵ * sᴮ);   # [i,]  

    H = [∂θ ∂β Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]
    I_meas = c .* (n * sᴮ)  # Measured current, ALBEDO added in later


    # ADD IN ALBEDO
    sᴵ_unscaled = sun_position(time) - pos;
    # ecl = eclipse_conical(-pos, sᴵ_unscaled) ####### NEED TO FIX TO +pos when updated
    # ecl = (ecl > 0.98) ? 1.0 : 0.0

    albedo_matrix, ignore = albedo(pos, sᴵ_unscaled, alb.refl)
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

# Pass in A by reference?
function prediction(x, w, dt, i, flag)
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
    β = x[5:7]; # Bias portion

    γ = w - β;     # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    θ = (nγ*dt);  
    r = γ/nγ;  # Make unit

    qp = qmult(q, [r*sin(θ/2); cos(θ/2)]); 
    
    skew = -hat(γ)

    R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues (for matrix exponential?)

    A = [R -dt*I(3); zeros(3,3) I(3)]; # Jacobian of f(x)

    A = [A          zeros(6, 3 * i); 
         zeros(3 * i, 6)     I(3*i)];

    xn = [qp; β]  # x at next step

    return xn, A
end
