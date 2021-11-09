####################################################################
#               DIODE CALIBRATION + MEKF                           #
####################################################################

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
    mag_calib_matrix # T
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

    new_state, new_covariance = mekf_sqrt(data.sat_state, sat.covariance, data.W, data.V,   # CHANGED data. -> sat.
                                        data.inertial_vecs, data.body_vecs, data.ang_vel, 
                                        data.current_meas, data.num_diodes, data.pos, data.dt, 
                                        data.time, data.albedo, data.mag_calib_matrix)

    # P = sat.covariance' * sat.covariance
    # state_alt, cov_alt = mekf(data.sat_state, P, data.W, data.V,   # CHANGED data. -> sat.
    #                                     data.inertial_vecs, data.body_vecs, data.ang_vel, 
    #                                     data.current_meas, data.num_diodes, data.pos, data.dt, 
    #                                     data.time, data.albedo)

    # if !(new_state ≈ state_alt) || !((new_covariance' * new_covariance) ≈ cov_alt)
    #     @infiltrate 
    # end

    # @btime mekf($data.sat_state, $sat.covariance, $data.W, $data.V,   
    #                 $data.inertial_vecs, $data.body_vecs, $data.ang_vel, 
    #                 $data.current_meas, $data.num_diodes, $data.pos, $data.dt, 
    #                 $data.time, $data.albedo)

    q, βgyro, βmag, c, α, ϵ = split_state(new_state, data.num_diodes)

    diodes = DIODES(c, α, ϵ) # Update with new estimates for calibration values, azimuth and elevation angles
    sat.diodes = diodes

    data.sat_state = new_state 

    sat.state = [q[:]; βgyro[:]; βmag[:]] # Only tracking the non-calibration states here
    sat.covariance = new_covariance

    return sat, data
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
                        system._num_diodes, # Number of photodiodes on the satellite 
                        generate_mag_calib_matrix(sat))   

    β₀_gyro = [0; 0; 0]
    β₀_mag  = [0; 0; 0]
    x₀ = [q; β₀_gyro; β₀_mag; sat.diodes.calib_values; sat.diodes.azi_angles; sat.diodes.elev_angles]

    if isnan(sat.covariance[1,1]) # Setting up covariance from scratch
        p = [σ_q * ones(3); σ_βgyro * ones(3); σ_βmag * ones(3); σ_c * ones(data.num_diodes); σ_α*ones(data.num_diodes); σ_ϵ*ones(data.num_diodes)].^2
        P₀ = diagm(p)
        P₀ = chol(P₀)
    else    # Setting up covariance after an eclipse
        p = [σ_q * ones(3); σ_βgyro * ones(3); σ_βmag * ones(3)].^2   # Reset Attitude and Bias covariance
        P_new = diagm(p)
        P₀ = sat.covariance # Should already be cholesky-ified
        P₀[1:9, 1:9] = chol(P_new)                    # Keep calibration covariances
    end

    # Not originally .^2 
    W = Diagonal([σ_orient * ones(3); σ_bias_gyro * ones(3); σ_bias_mag * ones(3); σ_cal * ones(data.num_diodes); σ_azi * ones(data.num_diodes); σ_ele * ones(data.num_diodes)]).^2
    V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)]).^2
    
    data.sat_state = x₀
    sat.covariance = P₀    
    data.W = chol(Matrix(W))
    data.V = chol(Matrix(V))

    return data
end

function mag_measurement(x, 𝐁ᴵ, i, T)
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
    q, βgyro, βmag, c, α, ϵ = split_state(x, i)

    ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion 
    𝐁ᴮ = ᴮQᴵ * 𝐁ᴵ;     # this is what the measurement would be given our estimated attitude

    B̌ᴮ = hat(𝐁ᴮ);    # Skew-symmetric matrix

    ∂θ = B̌ᴮ 
    ∂βgyro = zeros(3, 3)
    ∂βmag  = T
    ∂C = zeros(3, i)
    ∂α = zeros(3, i)
    ∂ϵ = zeros(3, i)    

    H = [∂θ ∂βgyro ∂βmag ∂C ∂α ∂ϵ]; # [3 x 6 + 3i]    #######
    y = 𝐁ᴮ[:]             # [3 x 1]

    return y, H
end

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

    q, βgyro, βmag, c, α, ϵ = split_state(x, i)

    ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (transposed to get I -> B)
    sᴮ = ᴮQᴵ * 𝐬ᴵ

    šᴮ= hat(sᴮ);  # Skew-symmetric form
    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    ndϵ = [(-sin.(ϵ).*cos.(α)) ((-sin.(ϵ).*sin.(α))) cos.(ϵ)]; # (With negative middle term, differing from the paper)
    
    ∂θ = (c .* n) * šᴮ;     # [i x 3]
    ∂βgyro = zeros(i, 3);       # [i x 3]  
    ∂βmag  = zeros(i, 3)
    ∂C = n * sᴮ;            # [i,]
    ∂α = c .* (ndα * sᴮ);   # [i,]
    ∂ϵ = c .* (ndϵ * sᴮ);   # [i,]  

    H = [∂θ ∂βgyro ∂βmag Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]
    I_meas = c .* (n * sᴮ)  # Measured current, ALBEDO added in later

    sᴵ_unscaled = sun_position(time) - pos;

    albedo_matrix, _ = albedo(pos, sᴵ_unscaled, alb.refl)

    for j = 1:i
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
    q, βgyro, βmag, c, α, ϵ = split_state(x, i)

    γ = w-βgyro;     # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    θ = (nγ*dt);  
    r = γ/nγ;  # Make unit

    qp = qmult(q, [r*sin(θ/2); cos(θ/2)]);   # Scalar-last quaternion prediciton
    
    skew = -hat(γ) 

    R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues (for matrix exponential?)

    # A = [R -dt*I(3); zeros(3,3) I(3)]; # Jacobian of non-calibration states  
    A = [R         -dt*I(3)   zeros(3,3);
         zeros(3,3) I(3)      zeros(3,3);
         zeros(3,3) zeros(3,3)      I(3)] 
    A = [A          zeros(9, 3 * i); 
        zeros(3 * i, 9)     I(3*i)];

    xn = [qp; βgyro; βmag; c; α; ϵ]; # x at next step

    return xn, A
end

function split_state(x, i)
    """ Helper function to split state into various components """
    x = x[:]
    q = x[1:4]
    βgyro = x[5:7]
    βmag  = x[8:10]
    c = x[11:(10+i)]  
    α = x[(11+i):(10+2*i)]
    ϵ = x[(11+2*i):end]

    return q, βgyro, βmag, c, α, ϵ
end
    
function mekf_sqrt(x, Pchol, cholW, cholV, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO, T) 

    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)                
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Prediction
    x_p, A = prediction(x, w, dt, _num_diodes); # State prediction
    Pchol_p = qrᵣ( [Pchol * A'; cholW] )

    # Measurement
    yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes, T)
    z_mag = 𝐁ᴮ - yp_mag 
    yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ,_num_diodes, pos, time, alb) 
    z_cur = y - yp_cur 

    C = [C_mag; C_cur]
    z = [z_mag[:]; z_cur[:]]

    Pyy_chol = qrᵣ( [Pchol_p * C'; cholV] )  
    L = ((Pchol_p' * Pchol_p * C') / Pyy_chol) / (Pyy_chol')
    # L₂ = (G \ ( (G' \ C) * F' * F) )'

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

    Pchol_next = qrᵣ( [Pchol_p*(I-L*C)'; cholV*L'] )

    return x_next, Pchol_next
end

function chol(M)
    # return cholesky(Symmetric(M)).U 
    return cholesky(Hermitian(M)).U
end

function qrᵣ(M)
    return qr(M).R 
end
