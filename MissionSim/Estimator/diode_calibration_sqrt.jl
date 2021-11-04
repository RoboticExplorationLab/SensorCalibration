####################################################################
#               DIODE CALIBRATION + MEKF                           #
####################################################################
using Infiltrator, BenchmarkTools, Test, MAT

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

    new_state, new_covariance = mekf_sqrt(data.sat_state, sat.covariance, data.W, data.V,   # CHANGED data. -> sat.
                                        data.inertial_vecs, data.body_vecs, data.ang_vel, 
                                        data.current_meas, data.num_diodes, data.pos, data.dt, 
                                        data.time, data.albedo)

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

    import sys
    sys.path.insert(0, "./") # Lets me import local files

    from numpy import arccos as acos # update

    from earthAlbedo import EarthAlbedo
    from .sqrtDiodeCalibrator import SqrtDiodeCalibrator
    import satelliteDynamics 
    import scipy.io

    _E_am0 = 1366.9


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
    σ_c = 0.15 # 0.2 
    σ_α = deg2rad(2.0) # 1.0 #2.0   # 2.0 is the σ used when generating these
    σ_ϵ = deg2rad(2.0) # 0.3 #1.0 

    if isnan(sat.covariance[1,1]) # Setting up covariance from scratch
        p = [σ_q * ones(3); σ_β * ones(3); σ_c * ones(data.num_diodes); σ_α*ones(data.num_diodes); σ_ϵ*ones(data.num_diodes)].^2
        P₀ = diagm(p)
        P₀ = chol(P₀)
    else    # Setting up covariance after an eclipse
        p = [σ_q * ones(3); σ_β * ones(3)].^2   # Reset Attitude and Bias covariance
        P_new = diagm(p)
        P₀ = sat.covariance # Should already be cholesky-ified
        P₀[1:6, 1:6] = chol(P_new)                    # Keep calibration covariances
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

    Q_diode = 1e-5   # Diode Noise  ##############################3
    σ_cal = Q_diode; σ_azi = Q_diode; σ_ele = Q_diode;
    σ_sunVec = deg2rad(3.0); σ_magVec = deg2rad(3.0); σ_curr = 0.05 #0.008; #3, 3, 0.005
    ##################################
    # NOTE (Not sure why σ_cal != σ_c... )
    # Not originally .^2 too
    W = Diagonal([σ_orient * ones(3); σ_bias * ones(3); σ_cal * ones(data.num_diodes); σ_azi * ones(data.num_diodes); σ_ele * ones(data.num_diodes)]).^2
    V = Diagonal([σ_magVec * ones(3); σ_curr * ones(data.num_diodes)]).^2
    
    data.sat_state = x₀
    sat.covariance = P₀    #chol(P₀) #Matrix(chol(P₀))
    data.W = chol(Matrix(W))
    data.V = chol(Matrix(V))

    return data
end

export mag_measurement
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
export current_measurement
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
    n = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
    ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
    ndϵ = [(-sin.(ϵ).*cos.(α)) ((-sin.(ϵ).*sin.(α))) cos.(ϵ)]; # (With negative middle term, differing from the paper)
    
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

# Make a "Rodrigues" helper function in the rotationFunctions.jl script?
export prediction
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

export split_state
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
    
export mekf_sqrt
function mekf_sqrt(x, Pchol, cholW, cholV, rᴵ, rᴮ, w, y, _num_diodes, pos, dt, time, alb::ALBEDO) 

    sᴵ = @view rᴵ[1,:];  𝐬ᴵ = sᴵ / norm(sᴵ)
    Bᴵ = @view rᴵ[2,:];  𝐁ᴵ = Bᴵ / norm(Bᴵ)                
    Bᴮ = @view rᴮ[2,:];  𝐁ᴮ = Bᴮ / norm(Bᴮ)

    # Prediction
    x_p, A = prediction(x, w, dt, _num_diodes); # State prediction
    # Pchol_p = qrᵣ( [Pchol * A'; chol(Matrix(W))] )  # Cholesky factor
    Pchol_p = qrᵣ( [Pchol * A'; cholW] )

    # Measurement
    yp_mag, C_mag = mag_measurement(x_p, 𝐁ᴵ, _num_diodes)
    z_mag = 𝐁ᴮ - yp_mag 
    yp_cur, C_cur = current_measurement(x_p, 𝐬ᴵ,_num_diodes, pos, time, alb) 
    z_cur = y - yp_cur 

    C = [C_mag; C_cur]
    z = [z_mag[:]; z_cur[:]]

    # Pyy_chol = qrᵣ( [Pchol_p * C'; chol(Matrix(V))] )  
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

    # temp = [Pchol_p*(I-L*C)'; chol(Matrix(V))*L']
    temp = [Pchol_p*(I-L*C)'; cholV*L']
    Pchol_next = qrᵣ( temp )

    return x_next, Pchol_next
end


# Helper functions 
export chol
function chol(M)
    # return cholesky(Symmetric(M)).U 
    return cholesky(Hermitian(M)).U
end

export qrᵣ
function qrᵣ(M)
    return qr(M).R 
end

# __init_diode_cal__()