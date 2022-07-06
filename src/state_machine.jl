# [src/state_machine.jl]

""" TO DO
 - Split by typeof(data)? 
 - Values for opt args
"""

# Given a state, generate measurements. Estimate/Control. Generate and return next step (maybe flip that to be first?)
#Consider splitting by typeof(data). Add in optional arguments.
""" step(sat_truth, sat_est, alb, x, t, dt, op_mode, flags, idx, progress_bar, T_orbit, data; use_albedo, initial_detumble_thresh, final_detumble_thresh, mag_ds_rate, 
            calib_cov_thres, mekf_cov_thres, σβ, σB, σ_gyro, σr, σ_current)

      Steps the system forward. Takes in a given next state, and generates measurements. The Estimator is used to update the 
    state estimate, and the controller generates a desired command when applicable. Integrates over a time step to get the next 
    time step, which is returned. 

      Time is updated, as well as flags. Whether or not there is an eclipse is predicted from the diode readings. 

    Arguments: (LOTS)
      - `sat_truth`: Satellite containing the true parameters for the satellite   |  SATELLITE
      - `sat_est`: Satellite containing best estimate of sat_truths' parameters   |  SATELLITE  
      - `alb`:  ALBEDO struct containing REFL and cell center info                |  ALBEDO
      - `x`:  Current state of the environment that is used to generate measurements  |  STATE
      - `t`:  Current time, as an Epoch                                           |  Epoch
      - `dt`: Time step                                                           |  Scalar 
      - `op_mode`: Current operation mode                                         |  Operation_mode 
      - `flags`: Struct used to track what mode to transition to                  |  FLAGS 
      - `idx`: How many time steps have elapsed since start                       |  Int 
      - `progress_bar`: Progress bar to display progress                          |  ProgressBar 
      - `T_orbit`: Time to complete a full orbit                                  |  Scalar 
      - `data`: Struct that contains whatever data is necessary for the Estimator |  (ANY)

      - `initial_detumble_thresh`: (Opt) Minimum angular velocity to terminate initial detumbling. Default is 15 deg/s
      - `final_detumble_thresh`:   (Opt) Minimum angular velocity to terminate final detumbling. Default is 8 deg/s
      - `mag_ds_rate`:  (Opt) Downsampling rate for the magnetometer calibration. Default is 60 
      - `calib_cov_thresh`:  (Opt) Threshold for the norm of the calibration terms in the covariance matrix to terminate diode_cal 
      - `mekf_cov_thresh`:  (Opt) Threshold for the norm of the non-calibration terms in the covariance matrix to terminate mekf
      (Others are noise parameters, see respective functions)
      
    Returns: 



"""
function step(sat_truth::SATELLITE, sat_est::SATELLITE, alb::ALBEDO, x::STATE{T}, t::Epoch, dt::Real, op_mode::Operation_mode, 
                flags::FLAGS, idx::Int, progress_bar, T_orbit, data; use_albedo = true, initial_detumble_thresh = deg2rad(15), final_detumble_thresh = deg2rad(8),  # was 25, 10 deg
                mag_ds_rate = 60, calib_cov_thres = 0.036, mekf_cov_thres = 0.004, σβ = 3.14e-5, σB = deg2rad(0.25), σ_gyro = 0.5e-4, 
                σr = 5e3, σ_current = 0.05 ) where {T}

    t += dt   # Update time

    ### Generate measurements

    truth, sensors, ecl, noise = generate_measurements(sat_truth, alb, x, t, dt; use_albedo = use_albedo, σB = σB, σ_gyro = σ_gyro, σr = σr, σ_current = σ_current);


    flags.magnetometer_calibrated && (sensors = correct_magnetometer(sat_est, sensors))  # Correct the magnetometer readings once it has been calibrated 

    ##### MAYBE ACCUMULATE AVERAGES DURING MAG CAL AND USE DEVIATION FROM MEAN?????            
    flags.in_sun = (flags.in_sun) ? norm(sensors.diodes ./ sat_est.diodes.calib_values) > 0.7 :
                                    norm(sensors.diodes ./ sat_est.diodes.calib_values) > 0.8 
    # flags.in_sun = ecl ### REMOVE !    

    ### Estimate & Control 
    next_mode = op_mode
    u = SVector{3, T}(zeros(3))
    notes = ""



    if op_mode == detumble      # -> mag_cal, mekf
        """ detumble -> mag_cal, mekf

            Starts by generating a control to slow the tumbling of the CubeSat. 
            There are two times that the sat is detumbling: the first is right after launch, 
            and it does an initial detumbling to allow for communication. The second occurs 
            after the photodiodes have been calibrated and a good estimate of the gyroscope 
            bias has been made, and allows for a more precise detumbling. 

            After the initial detumbling, we transition to calibrating the magnetometers. 
            After the final detumbling, we transition to the vanilla MEKF.

            No estimation occurs during this mode.
        """

        ctrl = DETUMBLER(sensors.gyro, sensors.magnetometer, dt)
        u    = generate_command(ctrl)

        # Gyro bias is estimated with diode calibration, so subtract it out 
        if flags.diodes_calibrated && (norm(sensors.gyro - sat_est.state.β) < final_detumble_thresh)
            flags.final_detumble = true

            if !(flags.in_sun)
                next_mode = chill
            else
                next_mode = mekf 
                data = MEKF_DATA()
                q = run_triad(sensors, sat_est, t, flags.in_sun; sᴮ_true = truth.ŝᴮ)  # x.q
                reset_cov!(sat_est; reset_calibration = false)
                sat_est = SATELLITE(sat_est.J, sat_est.magnetometer, sat_est.diodes, update_state(sat_est.state; q = q), sat_est.covariance)
            end

        elseif !flags.diodes_calibrated && (norm(sensors.gyro) < initial_detumble_thresh)
            flags.init_detumble  = true 

            if !(flags.in_sun)
                next_mode = chill
            else
                q = run_triad(sensors, sat_est, t, flags.in_sun; sᴮ_true = truth.ŝᴮ)  # x.q 
                sat_est = SATELLITE(sat_est.J, sat_est.magnetometer, sat_est.diodes, update_state(sat_est.state; q = q), sat_est.covariance)
                next_mode = mag_cal 
                # Bᴵ_pred = IGRF13(sensors.pos, t)
                # N_samples = Int(round(2 * T_orbit / (mag_ds_rate)))
                # data = MAG_CALIBRATOR(N_samples, vcat([sensors.magnetometer;]...), Bᴵ_pred);
                data = MEKF_DATA()
            end
        else 
            next_mode = detumble # Keep detumbling  
        end 

        notes = "Mode: $op_mode\t||̂ω||: $(norm(rad2deg.(sensors.gyro))) \t ||ω||: $(norm(rad2deg.(x.ω)))"


    elseif op_mode == mag_cal   # -> chill, diode_cal 
        """ mag_cal -> chill, diode_cal 

            Accumulates data over some 2 orbits and uses that data + Gauss-Newton 
            to estimate the magnetometer calibration parameters. To prevent unnecessarily
            large datasets, this downsamples at some prespecified rate, and then runs 
            the 'estimate' function when enough data has been gathered. 

            Transitions to 'diode_cal' if there is no eclipse; otherwise, it transitions 
            to 'chill' and waits.
        """ 

        # notes = "Mode: $op_mode \tSamples: $(data.idx[1])/$(data.N)"

        # Downsample 
        if true 
            if !(flags.in_sun)
                next_mode = chill 
            else 
                sat_est = Estimator.estimate(sat_est, sensors, data, alb, t, dt; use_albedo = use_albedo, calibrate = true)

                # Check if covariance of calibration states is low enough to fix
                if norm(sat_est.covariance[7:end, 7:end]) < 0.1 * calib_cov_thres 
                    next_mode = finished # detumble 
                    flags.diodes_calibrated = true 
                end

                notes = "mode: $op_mode \t||Σ|| = $(round(norm(sat_est.covariance[7:end, 7:end]), digits = 4)),\t ||Σs|| = $(round(norm(sat_est.covariance[7:9, 7:9]), digits = 4)) \t ||Σζ|| = $(round(norm(sat_est.covariance[10:12, 10:12]), digits = 4)) \t ||Σβ|| = $(round(norm(sat_est.covariance[13:15, 13:15]), digits = 4))"
            end
        end

    elseif op_mode == chill     # -> diode_cal, (mekf?)
        """ chill -> diode_cal 

            Temporary mode that is used as a waiting zone until something happens. 
            Right now, this is only called when diodes are being calibrated but 
            an eclipse is occuring. 

            Transitions to 'diode_cal' as soon as the eclipse is over. 
        """
        if flags.in_sun 
            next_mode = mag_cal # diode_cal 

            data = MEKF_DATA()
            q = run_triad(sensors, sat_est, t, flags.in_sun; sᴮ_true = truth.ŝᴮ)  # x.q 
            reset_cov!(sat_est; reset_calibration = true)
            sat_est = SATELLITE(sat_est.J, sat_est.magnetometer, sat_est.diodes, update_state(sat_est.state; q = q), sat_est.covariance)

        else 
            next_mode = chill
        end
        notes = "Mode: $op_mode\tEclipse: $ecl"



    elseif op_mode == diode_cal # -> chill, detumble (round 2)
        """ diode_cal -> chill, detumble (round 2)

            Runs the estimator for calibrating the diodes while estimating attitude and 
            gyroscope bias. Does not work during eclipses, so this checks and transitions 
            to 'chill' during eclipses.

            When the magnitude of the covariance of the calibration state C, α, and ϵ is 
            beneath some value, this transitions to 'detumble' for the final, more thorough 
            detumbling. 
        """
        if !(flags.in_sun) 
            next_mode = chill 
        else

            sat_est = Estimator.estimate(sat_est, sensors, data, alb, t, dt; use_albedo = use_albedo, calibrate_diodes = true)

            # Check if covariance of calibration states is low enough to fix
            if norm(sat_est.covariance[7:end, 7:end]) < calib_cov_thres 
                next_mode = detumble 
                flags.diodes_calibrated = true 
            end

            notes = "mode: $op_mode \t||Σ|| = $(norm(sat_est.covariance[7:end, 7:end])) \t ||ΣC|| = $(norm(sat_est.covariance[7:12, 7:12])) \t ||Σα|| = $(norm(sat_est.covariance[13:18, 13:18])) \t ||Σϵ|| = $(norm(sat_est.covariance[19:24, 19:24]))"
        end


    elseif op_mode == mekf      # -> finished
        """ mekf -> finished 

            Tracks the attitude and gyroscope bias of the CubeSat. Run once 
            all calibration is done, and once the covariance is small enough this 
            transitions to 'finished'

            Note that this must be preceeded with TRIAD to get the initial guess of q, as 
            well as reset_cov! and MEKF_DATA(), none of which are currently being done.
        """

        if false #!(flags.in_sun)
            next_mode = chill 
        else
            sat_est = Estimator.estimate(sat_est, sensors, data, alb, t, dt; use_albedo = use_albedo, calibrate_magnetometer = false)

            if norm(sat_est.covariance[1:6, 1:6]) < mekf_cov_thres 
                next_mode = finished 
            end 

            eq = norm((sat_truth.state.q ⊙ qconj(sat_est.state.q))[2:4])  # Quaternion error
            eβ = norm( sat_truth.state.β - sat_est.state.β)
            notes = "Mode: $op_mode \t||Σ|| = $(norm(sat_est.covariance[1:6, 1:6]))\t q Err: $eq\t β Err: $eβ"
        end

    elseif op_mode == finished 
        # Dont do anything
        notes = "Done!"
    else 
        @warn "Unrecognized mode!"
    end

    ### Update state 
    x⁺ = rk4(sat_truth.J, x, u, t, dt; σβ = σβ) # quat is normalized inside rk4 

    # Update sat_truth state to match 
    new_sat_state = SAT_STATE(x⁺.q, x⁺.β)
    sat_truth = SATELLITE(; J = sat_truth.J, mag = sat_truth.magnetometer, dio = sat_truth.diodes, sta = new_sat_state, cov = sat_truth.covariance)


    ### Display
    ProgressMeter.next!(progress_bar; showvalues = [(:Mode, op_mode), (:Iteration, idx), (:Notes, notes)])
    return sat_truth, sat_est, x⁺, t, next_mode, data, truth, sensors, ecl, noise
end



""" run_triad(sensors, sat_est, t, in_sun)

      Uses the current time and sensor measurements to estimate the 
    sun vector and magnetic field vector in body and inertial frames, and uses 
    TRIAD to initialize attitude estimate.

    Arguments:
      - `sensors`:  Sensor measurements at current time step    |  SENSORS
      - `sat_est`:  Current estimate of Satellite parameters    |  SATELLITE 
      - `t`:        Current time, as an Epoch                   |  Epoch 
      - `in_sun`:   (Opt) whether the satellite is eclipsed. Just for debugging. Default is true 

    Returns:
      - `q`:  Estimate of satellite attitude

"""
function run_triad(sensors::SENSORS{N, T}, sat_est::SATELLITE, t::Epoch, in_sun::Bool = true; sᴮ_true = nothing) where {N, T}
    (!in_sun) && @warn "run_triad should never be called if not in the sun!"

    sᴵ_est = sun_position(t) - sensors.pos     # Estimated sun vector 
    Bᴵ_est = IGRF13(sensors.pos, t)            # Estimated magnetic field vector
    ŝᴮ = estimate_sun_vector(sensors, sat_est.diodes; sᴮ_true = sᴮ_true)
    Bᴮ = sensors.magnetometer

    q, _ = triad(sᴵ_est, Bᴵ_est, ŝᴮ, Bᴮ)  # Write this function here too? 
    return SVector{4, T}(q) 
end

""" triad(r₁ᴵ, r₂ᴵ, r₁ᴮ, r₂ᴮ)

    Method for estimating the rotation matrix between two reference frames given one pair of vectors in each frame

    Arguments:
        - `r₁ᴵ`, `r₂ᴵ`: Pair of vectors in the Newtonian (inertial) frame     | [3,] each
        - `r₁ᴮ`, `r₂ᴮ`: Corresponding pair of vectors in body frame           | [3,]   each

    Returns:
        - `R`: A directed cosine matrix (DCM) representing the rotation     | [3 x 3]
                between the two frames 
        - `q`: A quaternion (scalar last) representing the rotation         | [4,]
                between the two frames
"""
function triad(r₁ᴵ,r₂ᴵ,r₁ᴮ,r₂ᴮ)

    𝐫₁ᴵ = r₁ᴵ / norm(r₁ᴵ)
    𝐫₂ᴵ = r₂ᴵ / norm(r₂ᴵ)
    𝐫₁ᴮ = r₁ᴮ / norm(r₁ᴮ)
    𝐫₂ᴮ = r₂ᴮ / norm(r₂ᴮ)

    t₁ᴵ = 𝐫₁ᴵ
    t₂ᴵ = cross(𝐫₁ᴵ, 𝐫₂ᴵ)/norm(cross(𝐫₁ᴵ,𝐫₂ᴵ));
    t₃ᴵ = cross(t₁ᴵ, t₂ᴵ)/norm(cross(t₁ᴵ,t₂ᴵ));

    Tᴵ = [t₁ᴵ[:] t₂ᴵ[:] t₃ᴵ[:]]

    t₁ᴮ = 𝐫₁ᴮ
    t₂ᴮ = cross(𝐫₁ᴮ, 𝐫₂ᴮ)/norm(cross(𝐫₁ᴮ,𝐫₂ᴮ));
    t₃ᴮ = cross(t₁ᴮ, t₂ᴮ)/norm(cross(t₁ᴮ,t₂ᴮ));

    Tᴮ = [t₁ᴮ[:] t₂ᴮ[:] t₃ᴮ[:]]

    R = Tᴵ * (Tᴮ')

    q = rot2quat(R);

    return q, R
end


# The Pseudo-inv method would probably be sketchy IRL because when the sun isn't illuminating 3+ diodes it would fail
# This method does not rely on knowledge of photodiode setup in advance
# Also, even with perfect estimates for the diode parameters this still gets ~2* error
function estimate_sun_vector(sens::SENSORS{N, T}, diodes::DIODES{N, T}; sᴮ_true = nothing) where {N, T}
    Ĩ = sens.diodes ./ diodes.calib_values 

    if false #size( Ĩ[abs.(Ĩ) .> 0.1], 1) < 3  # Less than three diodes are illuminated 
        @infiltrate
        return estimate_sun_vector2(sens, diodes)
    else
        return estimate_sun_vector_pinv(sens, diodes)
    end

end



function estimate_sun_vector_pinv(sens::SENSORS{N, T}, diodes::DIODES{N, T}) where {N, T}
    n(α, ϵ) = [sin(pi/2 - ϵ)*cos(α); sin(pi/2 - ϵ) * sin(α); cos(pi/2 - ϵ)]

    Ĩ    = sens.diodes ./ diodes.calib_values 
    azi  = diodes.azi_angles
    elev = diodes.elev_angles 
    
    ns = zeros(eltype(Ĩ),  6, 3) 
    for i = 1:6 
        ns[i, :] .= n(azi[i], elev[i])
    end 

    ŝ = (ns' * ns) \ (ns' * Ĩ)

    # REMOVE! 
    # if size(Ĩ[Ĩ .> 0.1], 1) < 3
    #     @infiltrate
    # end

    return ŝ / norm(ŝ)
end

function estimate_sun_vector2(sens::SENSORS{N, T}, diodes::DIODES{N, T}) where {N, T}

    sph2cart(α, ϵ, ρ) = [ρ * sin(pi/2 - ϵ)*cos(α); ρ * sin(pi/2 - ϵ) * sin(α); ρ * cos(pi/2 - ϵ)]
    
    sx, sy, sz = 0.0, 0.0, 0.0
    θ = deg2rad(-45)
    Ry = [cos(θ)  0  -sin(θ);
           0      1    0   ;
          sin(θ)  0   cos(θ)]

    _, _, z = sph2cart(diodes.azi_angles[1], diodes.elev_angles[1],  Ry * (sens.diodes[1] / diodes.calib_values[1]))    
    sz += z
    _, _, z = sph2cart(diodes.azi_angles[2], diodes.elev_angles[2],  Ry * (sens.diodes[2] / diodes.calib_values[2]))    
    sz += z

    _, y, _ = sph2cart(diodes.azi_angles[3], diodes.elev_angles[3],   sens.diodes[3] / diodes.calib_values[3])    
    sy += y
    _, y, _ = sph2cart(diodes.azi_angles[4], diodes.elev_angles[4],   sens.diodes[4] / diodes.calib_values[4])
    sy += y

    x, _, _ = sph2cart(diodes.azi_angles[5], diodes.elev_angles[5],  Ry * (sens.diodes[5] / diodes.calib_values[5]))    
    sx += x
    x, _, _ = sph2cart(diodes.azi_angles[6], diodes.elev_angles[6],  Ry * (sens.diodes[6] / diodes.calib_values[6]))    
    sx += x


    ŝᴮ = [sx, sy, sz]
    return SVector{3, T}(ŝᴮ / norm(ŝᴮ))
end

function estimate_sun_vector_old(sens::SENSORS{N, T}, sat_est::SATELLITE) where {N, T}
    """ Estimates a (unit) sun vector using the diode measurements 
            (Note that the equation is strange because a 45° rotation @ ŷ was used to avoid boundary problems with the elevation angles, bc azimuth isnt defined at ±Z) """

    n(α, ϵ) = [sin(pi/2 - ϵ)*cos(α);  sin(pi/2 - ϵ) * sin(α);  cos(pi/2 - ϵ)]
    
    # @warn "Estimate sun vector doesnt use the surface normals, or the updated angles!"
    if true # norm(sens.diodes) > eclipse_threshold  # If not eclipsed

        x₁ = (sens.diodes[1]/sat_est.diodes.calib_values[1])
        x₂ = (sens.diodes[2]/sat_est.diodes.calib_values[2])
        y₁ = (sens.diodes[3]/sat_est.diodes.calib_values[3])
        y₂ = (sens.diodes[4]/sat_est.diodes.calib_values[4])
        z₁ = (sens.diodes[5]/sat_est.diodes.calib_values[5])
        z₂ = (sens.diodes[6]/sat_est.diodes.calib_values[6])

        Is = [x₁, x₂, y₁, y₂, z₁, z₂]

        ns = zeros(6, 3)
        for i = 1:6
            # nᵢ = n(sat_est.diodes.azi_angles[i], sat_est.diodes.elev_angles[i])
            # ns[i, :] .= (Is[i] > 0.0) ? nᵢ  : zeros(3)
            ns[i, :] .= n(sat_est.diodes.azi_angles[i], sat_est.diodes.elev_angles[i])
        end
        sun_vec_est = (ns' * ns) \ (ns' * Is)

        return (sun_vec_est / norm(sun_vec_est))


        # sun_vec_est = [x₁ - x₂;
        #                 y₁ - y₂;
        #                 z₁ - z₂]
        # sun_vec_est /= norm(sun_vec_est)      

        # # Unrotate by 45* about +Y
        # θ = deg2rad(45)
        # Ry = [cos(θ)  0  -sin(θ);
        #        0      1    0   ;
        #       sin(θ)  0   cos(θ)]

        # sun_vec_est = Ry * sun_vec_est  

        # # Unrotate by 45* about +Y
        # sun_vec_est = [(x₁*cos(-pi/4) + z₁*cos(pi/4) + x₂*cos(3*pi/4) + z₂*cos(-3*pi/4));    
        #                 y₁ - y₂;
        #                (x₁*cos(3*pi/4) + z₁*cos(pi/4) + x₂*cos(-pi/4) + z₂*cos(-3*pi/4))] 

        # sun_vec_est /= norm(sun_vec_est)
    else
        # Assume that we are in eclipse  (this should never be called in eclipse though)
        sun_vec_est = [0; 0; 0]
    end
        
    return SVector{3, T}(sun_vec_est) # Unit - ŝᴮ
end


""" correct_magnetometer(sat, sensors)

      Uses the current estimate of the satellite magnetometer parameters to 
    correct the magnetometer reading.

    Arguments:
      - `sat`: Satellite with magnetometer parameters to use in calibrating (probably the sat_est)    |  SATELLITE 
      - `sensors`:  Set of sensor measurements at the current time step                               |  SENSORS 

    Returns:
      - `sensors`:  Set of sensor measurements with a corrected magnetometer reading                  |  SENSORS
"""
function correct_magnetometer(sat::SATELLITE, sensors::SENSORS{N, T}) where {N, T}

    a, b, c = sat.magnetometer.scale_factors 
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles 
    β = sat.magnetometer.bias

    T̂ = [a          0               0;
         b*sin(ρ)   b*cos(ρ)        0; 
         c*sin(λ)   c*sin(ϕ)*cos(λ) c*cos(ϕ)*cos(λ)  ]

    B̂ᴮ = SVector{3, T}(T̂ \ (sensors.magnetometer - β))

    
    return SENSORS(B̂ᴮ, sensors.diodes, sensors.gyro, sensors.pos)
end

""" Alternative form that uses the measured magnetic field vector in body frame rather than the sensors. """
function correct_magnetometer(sat::SATELLITE, B::SVector{3, T}) where {T}

    a, b, c = sat.magnetometer.scale_factors 
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles 
    β = sat.magnetometer.bias

    T̂ = [a          0               0;
         b*sin(ρ)   b*cos(ρ)        0; 
         c*sin(λ)   c*sin(ϕ)*cos(λ) c*cos(ϕ)*cos(λ)  ]

    B̂ᴮ = SVector{3, T}(T̂ \ (B - β))

    
    return B̂ᴮ
end

