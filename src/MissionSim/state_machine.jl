# Temporary, just used for debugging 
final_count = 0 
temp_count = 0

function update_operation_mode(flags::FLAGS, sens::SENSORS, system, albedo, current_data, t, satellite_estimate)


    # if flags.magnetometer_calibrated
    #     return finished, TRIVIAL(1.0), TRIVIAL(1.0), flags
    # end

    if flags.diodes_calibrated
        global temp_count += 1 
    else 
        global temp_count = 0
    end 
    if temp_count > 50
        return finished, TRIVIAL(1.0), TRIVIAL(1.0), flags 
    end
    
    # Prevent flickering due to noise
    if flags.in_sun 
        in_eclipse = norm(sens.diodes ./ satellite_estimate.diodes.calib_values) < 0.8
    else
        in_eclipse = norm(sens.diodes ./ satellite_estimate.diodes.calib_values) < 0.9
    end

    # Split sensor measurements
    ŝᴵ = sun_position(t) - sens.gps     # Estimated sun vector 
    B̂ᴵ = IGRF13(sens.gps, t)            # Estimated magnetic field vector
    Bᴮ = sens.magnetometer  
    
    𝐬̂ᴮ = in_eclipse ? SVector(0.0, 0.0, 0.0) : estimate_sun_vector(sens, satellite_estimate)
    ω = flags.diodes_calibrated ? (copy(sens.gyro) - satellite_estimate.state[5:7]) : copy(sens.gyro)

    if flags.magnetometer_calibrated
        Bᴮ = correct_mag_field(satellite_estimate, Bᴮ)  
    end


    ###########################################
    ###### FINISH CURRENT TASK FIRST ##########
    ###########################################

    if flags.detumbling 
        flags.detumbling = flags.diodes_calibrated ? (!check_if_finished(ω, deg2rad(10))) : (!check_if_finished(ω, deg2rad(25)))
        
        if flags.detumbling  # Still detumbling, continue
            mode, cont, est = detumble,  DETUMBLER(sens.gyro, Bᴮ, system._dt), TRIVIAL(0.0)
            return mode, cont, est, flags
        end

    elseif flags.calibrating
        if !flags.magnetometer_calibrated # Calibrating MAGNETOMETER
            if check_if_finished(satellite_estimate) 
                flags.calibrating, flags.magnetometer_calibrated = false, true 
                Bᴮ = correct_mag_field(satellite_estimate, sens.magnetometer)  # If it is now calibrated, correct Bᴮ 
            end       
        else  # Calibrating DIODES   
            if !in_eclipse # norm(sens.diodes) > eclipse_threshold # If still in sun
                if check_if_finished(satellite_estimate.covariance[7:end, 7:end], 0.07) # 01 0.007) #for non-sqrt
                    flags.calibrating, flags.diodes_calibrated = false, true
                end
            else
                flags.in_sun = false 
                flags.calibrating = false
            end
        end

        if flags.calibrating  # If still not done calibrating...
            if !flags.magnetometer_calibrated # Continue calibrating magnetometer
                mode, cont, est = mag_cal, TRIVIAL(0.0), MAG_CALIB(Bᴮ, B̂ᴵ)
            
            else 
                mode, cont = diode_cal, TRIVIAL(0.0)

                est = DIODE_CALIB(albedo,                       # Albedo
                                    current_data.sat_state,     # Satellite State
                                    # current_data.covariance,    # Filter Covariance 
                                    [ŝᴵ B̂ᴵ]',                   # Vectors in Inertial Frame
                                    [𝐬̂ᴮ Bᴮ]',                   # Vectors in Body Frame
                                    sens.gyro,                  # Angular Velocity, ω
                                    sens.diodes,                # Current Measurements, Ĩ
                                    sens.gps,                   # Position 
                                    current_data.W,             # 
                                    current_data.V,             #
                                    system._dt,                 # Time step, dt
                                    current_data.time,          #
                                    length(sens.diodes))        # Number of photodiodes

            end
            return mode, cont, est, flags, satellite_estimate
        end
    end


    ###########################################
    ###### SELECT NEXT TASK IF NEEDED #########
    ###########################################

    if ((flags.diodes_calibrated) && (norm(ω) > deg2rad(15.0))) || ((!flags.diodes_calibrated) && (norm(ω) > deg2rad(30.0)))
        mode = detumble
        cont = DETUMBLER(sens.gyro, Bᴮ, system._dt)  
        est = TRIVIAL(0.0)
        flags.detumbling = true

    elseif !flags.magnetometer_calibrated 
        mode = mag_cal
        cont = TRIVIAL(0.0)
        est = MAG_CALIB(Bᴮ, B̂ᴵ)
        flags.calibrating = true

    elseif !in_eclipse #norm(sens.diodes) > eclipse_threshold # (norm(𝐬̂ᴮ) > 0.9) # in sun -> Make consistent across everything!

        flags.in_sun = true 
        if !flags.diodes_calibrated
            ### CALIBRATE DIODES
            q, R = triad(ŝᴵ, B̂ᴵ, 𝐬̂ᴮ, Bᴮ) # Normalizes all input Vectors
            mode = diode_cal 
            cont = TRIVIAL(0.0)

            est = new_diode_calib(albedo, sens, system, q, satellite_estimate) # Init DIODE_CALIB
            est.inertial_vecs = [ŝᴵ B̂ᴵ]'
            est.body_vecs = [𝐬̂ᴮ Bᴮ]'
            est.current_meas = sens.diodes

            flags.calibrating = true
        else 
            mode = mekf 
            q, R = triad(ŝᴵ, B̂ᴵ, 𝐬̂ᴮ, Bᴮ)
            if isa(current_data, DIODE_CALIB)    # Finished diode cal during sunlit cycle 
                est = new_mekf_data(current_data)  
                # est = new_mekf_data(albedo, sens, system, q, satellite_estimate)
            elseif isa(current_data, MEKF_DATA)  # Continuing 
                est = current_data
            else    # Starting out of an eclipse 
                est = new_mekf_data(albedo, sens, system, q, satellite_estimate)
            end

            est.inertial_vecs = [ŝᴵ B̂ᴵ]'
            est.body_vecs = [𝐬̂ᴮ Bᴮ]'
            est.current_meas = sens.diodes
            est.pos = sens.gps
            est.ang_vel = sens.gyro

            cont = TRIVIAL(0.0)
            flags.calibrating = false
        end

    else 
        flags.in_sun = false 
        mode = chill
        cont = TRIVIAL(0.0)
        est = TRIVIAL(0.0)
        flags.calibrating = false
    end

    return mode, cont, est, flags, satellite_estimate
end



# TODO verify that this is still needed (also, cant I just do this with an R matrix? (YES) Also, does this work...?)
function estimate_sun_vector(sens::SENSORS, sat_est::SATELLITE)
    """ Estimates a (unit) sun vector using the diode measurements 
            (Note that the equation is strange because a 45° rotation @ ŷ was used to avoid boundary problems with the elevation angles) """

    if true # norm(sens.diodes) > eclipse_threshold  # If not eclipsed
        x₁ = (sens.diodes[1]/sat_est.diodes.calib_values[1])
        x₂ = (sens.diodes[2]/sat_est.diodes.calib_values[2])
        y₁ = (sens.diodes[3]/sat_est.diodes.calib_values[3])
        y₂ = (sens.diodes[4]/sat_est.diodes.calib_values[4])
        z₁ = (sens.diodes[5]/sat_est.diodes.calib_values[5])
        z₂ = (sens.diodes[6]/sat_est.diodes.calib_values[6])

        # sun_vec_est = [x₁ - x₂;
        #                 y₁ - y₂;
        #                 z₁ - z₂]

        sun_vec_est = [(x₁*cos(-pi/4) + z₁*cos(pi/4) + x₂*cos(3*pi/4) + z₂*cos(-3*pi/4));    
                        y₁ - y₂;
                       (x₁*cos(3*pi/4) + z₁*cos(pi/4) + x₂*cos(-pi/4) + z₂*cos(-3*pi/4))] 

        sun_vec_est /= norm(sun_vec_est)
    else
        # Assume that we are in eclipse  (this should never be called in eclipse though)
        sun_vec_est = [0; 0; 0]
    end
        
    return sun_vec_est # Unit - ŝᴮ
end

function correct_mag_field(sat, B̃ᴮ)
    """ Uses magnetometer calibration stuff to fix Bᴮ using sat ESTIMATE! """
    # Generate T, Bias from sat 
    a, b, c = sat.magnetometer.scale_factors 
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles 
    β = sat.magnetometer.bias

    T̂ = [a          0               0;
         b*sin(ρ)   b*cos(ρ)        0; 
         c*sin(λ)   c*sin(ϕ)*cos(λ) c*cos(ϕ)*cos(λ)  ]

    B̂ᴮ = T̂^(-1) * (B̃ᴮ - β)

    return B̂ᴮ # Corrected mag field in body frame -> not unit
end

function triad(r₁ᴵ,r₂ᴵ,r₁ᴮ,r₂ᴮ)
    """
        Method for estimating the rotation matrix between two reference frames given one pair of vectors in each frame

        Arguments:
            - r₁ᴵ, r₂ᴵ: Pair of vectors in the Newtonian (inertial) frame     | [3,] each
            - r₁ᴮ, r₂ᴮ: Corresponding pair of vectors in body frame           | [3,]   each

        Returns:
            - R: A directed cosine matrix (DCM) representing the rotation     | [3 x 3]
                    between the two frames 
            - q: A quaternion (scalar last) representing the rotation         | [4,]
                    between the two frames
    """

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

    q = q_from_DCM(R);

    return q, R
end


# For detumbler and diodes
function check_if_finished(measurement, thresh)
    """ Compares the norm of the measurement to a provided threshold to determine if it is finished. 
            Used primarily for detumbler (ang vel) and diode calibration (covariance) """
    return (norm(measurement) < thresh)
end

# For magnetometer 
function check_if_finished(sat_est::SATELLITE)
    """ Checks if the magnetometer has run. If it has, verifies that it did not result in any major non-orthogonality angles, which invalidate the solution """
    if !check_if_run()
        return false
    end

    ρ_est, λ_est, ϕ_est = sat_est.magnetometer.non_ortho_angles 
    if (abs(ρ_est) > pi/3) || (abs(λ_est) > pi/3) || (abs(ϕ_est) > pi/3)
        println("Error with major non-ortho angles!")
        return false
    else
        println("Finished MAGNETOMETER")
        return true
    end
end
