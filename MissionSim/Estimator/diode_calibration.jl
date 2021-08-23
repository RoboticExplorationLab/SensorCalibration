####################################################################
#               DIODE CALIBRATION + MEKF                           #
####################################################################
using Infiltrator, BenchmarkTools

mutable struct DIODE_CALIB 
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

function estimate_vals(sat::SATELLITE, data::DIODE_CALIB)
    """
        Calls a single iteration of a multiplicative extended Kalman filter 
                and updates satellite estimates, including calibration parameters

        Arguments:
            - sat:  Struct containing current estimate of satellite values      |  SATELLITE 
            - data: Struct containing necessary information to run MEKF         |  -----

        Returns:
            - sat:  Updated satellite struct containing new estimates           |  SATELLITE
            - data: Updated struct containing necessary MEKF information        |  -----
    """  
    data.time += data.dt   
    new_state, new_covariance = mekf(data.sat_state, sat.covariance, data.W, data.V,   # CHANGED data. -> sat.
                                        data.inertial_vecs, data.body_vecs, data.ang_vel, 
                                        data.current_meas, data.num_diodes, data.pos, data.dt, 
                                        data.time, data.albedo)

    # @btime mekf($data.sat_state, $sat.covariance, $data.W, $data.V,   
    #                 $data.inertial_vecs, $data.body_vecs, $data.ang_vel, 
    #                 $data.current_meas, $data.num_diodes, $data.pos, $data.dt, 
    #                 $data.time, $data.albedo)

    q, β, c, α, ϵ = split_state(new_state, data.num_diodes)

    diodes = DIODES(c, α, ϵ) # Update with new estimates for calibration values, azimuth and elevation angles
    sat.diodes = diodes

    data.sat_state = new_state 

    sat.state = [q[:]; β[:]] # Only tracking the non-calibration states here
    sat.covariance = new_covariance

    return sat, data 
end

function __init_diode_cal__()
    py"""
    import numpy as np 
    # import ulab.numpy as np 
    _E_am0 = 1366.9

    def diag(v): # Copy of numpy, because ulab doesnt have it, but also abbreviated (?)
        n = v.shape[0]
        result = np.zeros((n,n), v.dtype)
        result[:n].flat[0::n+1] = v 
        return result

    def split_state(x, i):
        q = x[0:4]
        bias = x[4:7]
        c = x[7:(7+i)]
        azi = x[(7+i):(7+2*i)]
        elev = x[(7+2*i):]

        return q, bias, c, azi, elev
    
    def mag_measurement(x, B_inert, i):
        q, b, c, a, e = split_state(x, i)

        bQi = dcm_from_q(q).T 
        B_body = bQi @ B_inert 

        dtheta = hat(B_body)
        dbias  = np.zeros([3,3])
        dc     = np.zeros([3,i])
        da     = np.zeros([3,i])
        de     = np.zeros([3,i])

        H = np.hstack([dtheta, dbias, dc, da, de])
        y = B_body 

        return y, H 

    def prediction(x, w, dt, i):
        q, b, c, a, e = split_state(x, i)
        
        gamma = w - b   # Adjusted angular velocity (measured - bias)
        gamma_mag = np.linalg.norm(gamma)
        gamma_skew = -hat(gamma) / gamma_mag

        theta = gamma_mag * dt
        r = gamma / gamma_mag  # unit 

        quat_pred = qmult(q, np.hstack([r*np.sin(theta/2), np.cos(theta/2)]))

        R = (np.eye(3) + (gamma_skew * np.sin(theta)) + (gamma_skew @ gamma_skew)*(1 - np.cos(theta)))

        A = np.vstack([ np.hstack([R, -dt*np.eye(3)]),
                        np.hstack([np.zeros([3,3]), np.eye(3)]) ])

        A = np.vstack([ np.hstack([A, np.zeros([6, 3*i])]),
                        np.hstack([np.zeros([3*i, 6]), np.eye(3*i)])   ])

        x_pred = np.hstack([quat_pred, b, c, a, e])

        return x_pred, A
    
    def current_measurement(x, s_inert, i, pos, time, alb_refl, alb_cell_centers_ecef):
        q, b, c, a, e = split_state(x, i)
        bQi = dcm_from_q(q).T 
        s_body = bQi @ s_inert 

        s_body_skew = hat(s_body)
        surf_norm = np.array([np.cos(e)*np.cos(a),
                              np.cos(e)*np.sin(a),
                              np.sin(e)])

        dnda = np.array([-np.cos(e)*np.sin(a),
                          np.cos(e)*np.cos(a),
                          np.zeros(size(a))])
        dnde = np.array([-np.sin(e)*np.cos(a),
                         -np.sin(e)*np.sin(a), # negated
                          np.cos(e)])

        dtheta = (c * surf_norm) * s_body_skew 
        dbias  = np.zeros([i, 3])
        dc     = surf_norm @ s_body
        da     = c * (dnda @ s_body)
        de     = c * (dnde @ s_body)

        H = np.hstack([dtheta, dbias, diag(dc), diag(da), diag(de)])
        I = c * (surf_norm @ s_body)

        s_inert_unscaled = sun_position(time) - pos

        albedo_matrix, _ = albedo(pos, s_inert_unscaled, alb_refl)

        diode_albedos = compute_diode_albedo(albedo_matrix, alb_cell_centers_ecef, surf_norm, pos)
        diode_albedos = c * diode_albedos / _E_am0
        I = I + diode_albedos 


        I[I < 0] = 0 
        H[I < 0] = 0
        y = I 
        return y, H

    
    def sun_position(time): # From SatelliteDynamics? Make these into a .py for others 
        pass

    def compute_diode_albedo(alb_mat, cell_centers, surf_norm, pos):
        pass

    def albedo(pos, si, refl):
        pass

    """
end

function initialize(albedo, state, system) 
    """
        Temporary function for testing. Used to initialize a struct with true values rather than estimates
    """
        @info "Using true state to initialize diode calibration!"
        sat_state = [state[7:10][:]; zeros(3 + 3*system._num_diodes)[:]] 
        d = DIODE_CALIB(albedo,         # ALBEDO 
                        sat_state,      # Satellite state ([(q⃗, q₀) β C α ϵ], NOT the same as state in simulation)
                        # 0.0, # Just empty stuff to be filled in 
                        0.0,            # Ref vectors in inertial frame (rᴵ)
                        0.0,            # Ref vectors in body frame (rᴮ)
                        state[11:13],   # Angular velocity  
                        0.0,            # Photodiode current measurements
                        state[1:3],     # Position
                        0.0,            # Update with W 
                        0.0,            # Update with V 
                        system._dt,     # Time step
                        system._epc,    # Initial time
                        system._num_diodes)   # Number of photodiodes on the satellite 

        return d
end

function new_diode_calib(albedo, sens::SENSORS, system, q, sat) 
    """
        Generates a new struct containing necessary information for the diode calibration MEKF

        Arguments:
            - albedo: ALBEDO struct containing refl and cell centers      |  ALBEDO
            - sens: SENSORS struct containing sensor measurements         |  SENSORS 
            - system: Set of system constants needed for Simulation       |  List
            - q: Current estimate of satellite attitude                   |  [4,] ([q⃗; q₀] -> Scalar last)
            - sat: Struct containing current satellite estimate           |  SATELLITE

        Returns:
            - data: DIODE_CALIB struct used to run the calibration MEKF   |  DIODE_CALIB
    """     
    sat_state = zeros(7 + 3 * system._num_diodes)
    data = DIODE_CALIB(albedo,        # ALBEDO 
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

    β₀ = [0; 0; 0]
    x₀ = [q; β₀; sat.diodes.calib_values; sat.diodes.azi_angles; sat.diodes.elev_angles]

    σ_q = (10*pi/180) 
    σ_β = (10*pi/180)
    σ_c = 0.2 #0.15
    σ_α = 1.0 #deg2rad(2.0) #1.0;
    σ_ϵ = 0.3 #deg2rad(2.0) #0.3

    if isnan(sat.covariance[1,1]) # Setting up covariance from scratch
        p = [σ_q * ones(3); σ_β * ones(3); σ_c * ones(data.num_diodes); σ_α*ones(data.num_diodes); σ_ϵ*ones(data.num_diodes)].^2
        P₀ = diagm(p)
    else    # Setting up covariance after an eclipse
        p = [σ_q * ones(3); σ_β * ones(3)].^2   # Reset Attitude and Bias covariance
        P_new = diagm(p)
        P₀ = sat.covariance 
        P₀[1:6, 1:6] = P_new                    # Keep calibration covariances
    end


    ##### VERIFY THIS SECTION ########
    estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
                        gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
                        velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
                        accel_bias_instability = 6)      # in microG

    Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  
    σ_orient = sqrt(Q_gyro);

    Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)   # This is super small
    σ_bias = sqrt(Q_bias)

    Q_diode = 1e-6   # Diode Noise  ##############################3
    σ_cal = Q_diode; σ_azi = Q_diode; σ_ele = Q_diode;
    σ_sunVec = deg2rad(5.0); σ_magVec = deg2rad(5.0); σ_curr = 0.008; #3, 3, 0.005
    ##################################

    W = Diagonal([σ_orient * ones(3); σ_bias * ones(3); σ_cal * ones(data.num_diodes); σ_azi * ones(data.num_diodes); σ_ele * ones(data.num_diodes)])
    V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)])
    
    data.sat_state = x₀
    sat.covariance = P₀
    data.W = W 
    data.V = V

    return data
end

function mekf(x, P, W, V, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO)
    """
        Runs a single step of a multiplicative extended Kalman filter

        Arguments:
            - x: Current state + calibration estimates [(q⃗, q₀) β C α ϵ]  |  [7 + 3i,]
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
    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)                    
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Prediction
    x_p, A = prediction(x, w, dt, _num_diodes); # State prediction

    # ##### REMOVE #####
    # x_p_alt, A_alt = py"prediction"(x, w, dt, _num_diodes)
    # if (norm(x_p - x_p_alt) > 0.01) || (norm(A - A_alt) > 0.01)
    #     @infiltrate 
    # end
    # ##################

    P_p = A*P*A' + W;                           # Covariance prediction 

    # Measurement
    yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes)

    # ##### REMOVE ##### 
    # yp_mag_alt, C_mag_alt = py"mag_measurement"(x_p, 𝐁ᴵ, _num_diodes)
    # if (norm(yp_mag - yp_mag_alt) > 0.01) || (norm(C_mag - C_mag_alt) > 0.01)
    #     @infiltrate 
    # end
    # ##################

    z_mag = 𝐁ᴮ - yp_mag 
    yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ,_num_diodes, pos, time, alb) 
    z_cur = y - yp_cur 

    C = [C_mag; C_cur]
    z = [z_mag[:]; z_cur[:]]

    # Innovation   
    Vk = Diagonal(V) 
    S = C*P_p*C' + Vk; 

    # Kalman Gain
    L = P_p * C' * S^(-1); 

    # Update
    dx = L*z;          
    dϕ   =  @view dx[1:3]; 
    drest = @view dx[4:end]

    θ = (norm(dϕ));
    rTemp = dϕ / θ; 
    
    dq = [rTemp*sin(θ/2); cos(θ/2)];

    x_next = deepcopy(x)
    x_next[1:4] = qmult(view(x_p, 1:4), dq); 
    x_next[5:end] = view(x_p, 5:length(x_p)) + drest;
    
    T = (I - L*C) 
    P_next = T * P_p * T' + L*Vk*L';  
    # P_next = (I(size(P,1)) - L*C) * P_p * (I(size(P,1)) - L*C)' + L*Vk*L';  

    return x_next, P_next
end

function mag_measurement(x, 𝐁ᴵ, i)
    """
        Generates what we would expect our measured magnetic field vector would be in the body frame 
            given our current estimated attitude

        Arguments:
            - x:  Current state of the satellite [(q⃗, q₀) β C α ϵ]             | [7 + 3i,] 
            - 𝐁ᴵ: Unit mag field vector in the newtonian (inertial) frame      |  [3,]
            - i:  Number of photodiodes                                        |  Int

        Returns:
            - y: Unit mag field vector estimate corresponding to state         |  [3,]
            - H: Jacobian of y with respect to x                               |  [3 x 6 + 3i]
                        dy/dx = [dy/dϕ; dy/dβ; ...]                             
    """
    q, β, c, α, ϵ = split_state(x, i)

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

# NO LONGER ESTIMATES ECL!
function current_measurement(x, 𝐬ᴵ, i, pos, time, alb::ALBEDO)  
    """
        Generates what we would expect our measured current values to be given our estimate of state

        Arguments:
            - x:  Current state of the satellite [(q⃗, q₀) β C α ϵ]             | [7 + 3i,] 
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
    q, β, c, α, ϵ = split_state(x, i)

    ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (transposed to get I -> B)
    sᴮ = ᴮQᴵ * 𝐬ᴵ

    šᴮ= hat(sᴮ);  # Skew-symmetric form
    n   = [ cos.(ϵ).*cos.(α)    cos.(ϵ).*sin.(α)      sin.(ϵ)];  # [i x 3]
    ndα = [(-cos.(ϵ).*sin.(α))  (cos.(ϵ).*cos.(α))    zeros(size(α))];
    ndϵ = [(-sin.(ϵ).*cos.(α))  ((-sin.(ϵ).*sin.(α))) cos.(ϵ)]; # (With negative middle term, differing from the paper)
    
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
        # surface_normal = [cos(ϵ[j])*cos(α[j]) cos(ϵ[j])*sin(α[j]) sin(ϵ[j])]     # Photodiode surface normal
        
        diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, view(n, j, :), pos)
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

function prediction(x, w, dt, i)
    """
        Predicts the next state and covariance using current state, angular velocity, and time step
            (essentially assumes a small rotation in attitude and that the other states are constant)

        Arguments:
            - x:  Current state of the satellite [(q⃗, q₀) β C α ϵ]              | [7 + 3i,] 
            - w:  Current angular velocity estimate of the satellite            | [3,]
            - dt: Time step of the Simulation                                   | Scalar
            - i:  Number of photodiodes                                         |  Int

        Returns:
            - xn: Predicted next state  [(q⃗, q₀) β C α ϵ]                       | [7 + 3i,] 
            - H:  Jacobian of state x with respect to itself                    |  [i x 6 + 3i]
                    (Note that the quaternions have been replaced with 3 param)
                        dx/dx = [dϕ/dϕ; dϕ/dβ; ...]                             
    """
    q, β, c, α, ϵ = split_state(x, i)

    γ = w-β;     # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    θ = (nγ*dt);  
    r = γ/nγ;  # Make unit

    qp = qmult(q, [r*sin(θ/2); cos(θ/2)]);   # Scalar-last quaternion prediciton
    
    skew = -hat(γ) 

    R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues (for matrix exponential?)

    A = [R -dt*I(3); zeros(3,3) I(3)]; # Jacobian of non-calibration states   
    A = [A          zeros(6, 3 * i); 
        zeros(3 * i, 6)     I(3*i)];

    xn = [qp; β; c; α; ϵ]; # x at next step

    return xn, A
end

function split_state(x, i)
    """ Helper function to split state into various components """
    x = x[:]
    q = x[1:4]
    β = x[5:7]
    c = x[8:(7+i)]  
    α = x[(8+i):(7+2*i)]
    ϵ = x[(8+2*i):end]

    return q, β, c, α, ϵ
end
    



##### TESTING

function mekf_iterative(x, P, W, V, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO)
    """
        Runs a single step of a multiplicative extended Kalman filter

        Arguments:
            - x: Current state + calibration estimates [(q⃗, q₀) β C α ϵ]  |  [7 + 3i,]
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
    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)
    𝐬ᴮ = @view rᴮ[1,:];  𝐬ᴮ = 𝐬ᴮ / norm(𝐬ᴮ)                    
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Prediction
    x_p, A = prediction(x, w, dt, _num_diodes); # State prediction
    P_p = A*P*A' + W;                           # Covariance prediction 

    # Measurement
    yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes)
    yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ,_num_diodes, pos, time, alb) 
    C = [C_mag; C_cur]
    yyp = [yp_mag[:]; yp_cur[:]]
    yy = [𝐁ᴮ[:]; y[:]]
    Vk = Diagonal(V) 

    x_p_og = copy(x_p)
    z_hat = yyp

    for iter = 1:3 
        yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes)
        yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ,_num_diodes, pos, time, alb) 
        C = [C_mag; C_cur]

        S = C*P_p*C' + Vk # P\_yy
        L = P_p * C' * S^(-1); # K
        dx = L*(yy - yyp);          
        dϕ = dx[1:3]; 
        drest = dx[4:end]

        θ = (norm(dϕ));
        rTemp = dϕ / θ; 
    
        dq = [rTemp*sin(θ/2); cos(θ/2)];

        x_next = deepcopy(x)
        x_next[1:4] = qmult(x_p_og[1:4], dq); 
        x_next[5:end] = x_p_og[5:end] + drest;
        x_p = x_next
    end
    # Innovation   
    S = C*P_p*C' + Vk;  # P\_yu
    # Kalman Gain
    L = P_p * C' * S^(-1); 
    z = yy - yyp


    # Update
    dx = L*z;          
    dϕ = dx[1:3]; 
    drest = dx[4:end]

    θ = (norm(dϕ));
    rTemp = dϕ / θ; 
    
    dq = [rTemp*sin(θ/2); cos(θ/2)];

    x_next = deepcopy(x)
    x_next[1:4] = qmult(x_p[1:4], dq); 
    x_next[5:end] = x_p[5:end] + drest;
    
    P_next = (I(size(P,1)) - L*C) * P_p * (I(size(P,1)) - L*C)' + L*Vk*L';  

    return x_next, P_next
end

# May need to be adapted for the multiplicative portion
function mekf_sequential(x, P, W, V, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO)
    """
        Runs a single step of a multiplicative extended Kalman filter

        Arguments:
            - x: Current state + calibration estimates [(q⃗, q₀) β C α ϵ]  |  [7 + 3i,]
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
    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)
    𝐬ᴮ = @view rᴮ[1,:];  𝐬ᴮ = 𝐬ᴮ / norm(𝐬ᴮ)                    
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Prediction
    x̂ₖ⁻, A = prediction(x, w, dt, _num_diodes); # State prediction
    Pₖ⁻ = A*P*A' + W;                           # Covariance prediction 

    
    # Set 1 - Magnetometer
        # Measurement
        zₐ = 𝐁ᴮ
        ẑₐ, Hₐ = mag_measurement(x̂ₖ⁻, 𝐁ᴵ, _num_diodes)
  
        # Innovation   
        Rₐ = @view V[1:3, 1:3] 
        S = Hₐ*Pₖ⁻*Hₐ' + Rₐ; 

        Kₐ = Pₖ⁻ * Hₐ' * S^(-1);  # Kalman Gain

        # Update
        dx = Kₐ*(zₐ - ẑₐ);    

        dϕ = dx[1:3]; 
        drest = dx[4:end]

        θ = (norm(dϕ));
        rTemp = dϕ / θ; 
        dq = [rTemp*sin(θ/2); cos(θ/2)];

        x̂ₖᵃ = deepcopy(x)
        x̂ₖᵃ[1:4] = qmult(x̂ₖ⁻[1:4], dq); 
        x̂ₖᵃ[5:end] = x̂ₖ⁻[5:end] + drest;
        
        Iₚ = I(size(P,1))
        Pₖᵃ = (Iₚ - Kₐ*Hₐ) * Pₖ⁻ * (Iₚ - Kₐ*Hₐ)' + Kₐ*Rₐ*Kₐ';  

    # Set 2 - Diodes
        zᵦ = y
        ẑᵦ, Hᵦ = current_measurement(x̂ₖᵃ, 𝐬ᴵ,_num_diodes, pos, time, alb) 
        Rᵦ = @view V[4:end, 4:end] 
        S = Hᵦ * Pₖᵃ * Hᵦ' + Rᵦ 
        Kᵦ = Pₖᵃ * Hᵦ' * S^(-1)

        # Update
        dx = Kᵦ*(zᵦ - ẑᵦ);          
        dϕ = dx[1:3]; 
        drest = dx[4:end]

        θ = (norm(dϕ));
        rTemp = dϕ / θ; 
        dq = [rTemp*sin(θ/2); cos(θ/2)];

        x_next = deepcopy(x)
        x_next[1:4] = qmult(x̂ₖ⁻[1:4], dq); 
        x_next[5:end] = x̂ₖ⁻[5:end] + drest;
        
        P_next = (Iₚ - Kᵦ*Hᵦ) * Pₖᵃ * (Iₚ - Kᵦ*Hᵦ)' + Kᵦ*Rᵦ*Kᵦ'; 

    return x_next, P_next
end

__init_diode_cal__()