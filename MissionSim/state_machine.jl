

# Should albedo be in system/params or something?
temp_count = 0
final_count = 0 

# Should we prevent it from getting stuck in "detumble" endlessly?

println("Using sᴮ truth")
function update_operation_mode(flags::FLAGS, sens::SENSORS, system, albedo, current_data, t, satellite_estimate, sᴮ_truth)
    #flags:= {in_sun, magnetometer_calibrated, diodes_calibrated, calibrating, detumbling}
    # Don't like data_structure names for estimators, or order of DIODE_CALIB

    ŝᴵ = sun_position(t) - sens.gps     # Estimated sun vector 
    B̂ᴵ = IGRF13(sens.gps, t)            # Estimated magnetic field vector
    𝐬̂ᴮ = estimate_sun_vector(sens, satellite_estimate)
    Bᴮ = sens.magnetometer      # SHOULD BE A tilde

    if flags.magnetometer_calibrated
        Bᴮ = correct_mag_field(satellite_estimate, Bᴮ)  # Should be a hat?
        # return finished, TRIVIAL(1.0), TRIVIAL(1.0), flags
    end

    if flags.diodes_calibrated && flags.magnetometer_calibrated
        global final_count += 1
        if final_count > 5000
            return finished, TRIVIAL(1.0), TRIVIAL(1.0), flags
        end
    end




    ###########################################
    ###### FINISH CURRENT TASK FIRST ##########
    ###########################################

    if flags.detumbling 
        if flags.diodes_calibrated
            temp = deepcopy(sens)
            temp.gyro = temp.gyro - satellite_estimate.state[5:7] # current_data.sat_state[5:7]
            flags.detumbling = !check_if_finished(sens, 10)  # detumbler - 5 °/sec, is detumbling if *not* finished
        else
            flags.detumbling = !check_if_finished(sens, 25) # detumbler - 
        end
        
        if flags.detumbling  # Still detumbling, continue
            mode, cont, est = detumble,  DETUMBLER(sens.gyro, Bᴮ, system._dt), TRIVIAL(0.0)
            return mode, cont, est, flags
        end

    elseif flags.calibrating
        if !flags.magnetometer_calibrated # Calibrating MAGNETOMETER
            if check_if_finished(satellite_estimate) 
                flags.calibrating, flags.magnetometer_calibrated = false, true 
                Bᴮ = correct_mag_field(satellite_estimate, sens.magnetometer)  # If it is now calibrated, correct Bᴮ 
                # return finished, TRIVIAL(1.0), TRIVIAL(1.0), flags 
            end       
        else #if !flags.diodes_calibrated # Calibrating DIODES   
            # DONT LIKE THIS METHOD FOR CHECKING FINISHED!
            if sum(abs.(sens.diodes)) > eclipse_threshold # If still in sun
                temp = 0
                if check_if_finished(current_data, satellite_estimate) # diodes -> Do we need to enforce one sunlight-orbit...?
                    flags.calibrating, flags.diodes_calibrated = false, true
                end
            else
                flags.in_sun = false 
                flags.calibrating = false
            end
        # else
        #     tttt = 5
        end

        if flags.calibrating  # If still not done calibrating...
            if !flags.magnetometer_calibrated # Continue calibrating magnetometer
                mode, cont, est = mag_cal, TRIVIAL(0.0), MAG_CALIB(Bᴮ, B̂ᴵ)
            
            else #if !flags.diodes_calibrated # Continue calibrating diodes
                mode, cont = diode_cal, TRIVIAL(0.0)

                est = DIODE_CALIB(albedo,                       # Albedo
                                    current_data.sat_state,     # Satellite State
                                    # current_data.covariance,    # Filter Covariance 
                                    [ŝᴵ B̂ᴵ]',                   # Vectors in Inertial Frame
                                    sens.gyro,                  # Angular Velocity, ω
                                    [𝐬̂ᴮ Bᴮ]',                   # Vectors in Body Frame
                                    sens.diodes,                # Current Measurements, Ĩ
                                    current_data.W,             # 
                                    current_data.V,             #
                                    system._dt,                 # Time step, dt
                                    current_data.time,          #
                                    length(sens.diodes),        # Number of photodiodes
                                    sens.gps,                   # Position 
                                    false)                      # First pass 
            end
            return mode, cont, est, flags, satellite_estimate
        end
    end


    ###########################################
    ###### SELECT NEXT TASK IF NEEDED #########
    ###########################################
    
    if flags.diodes_calibrated
        ω = sens.gyro - satellite_estimate.state[5:7] 
    else
        ω = sens.gyro
    end    

    sᴵ_unscaled = sun_position(t) - sens.gps;
    ecl = eclipse_conical(-sens.gps, sᴵ_unscaled)


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

    elseif (norm(𝐬̂ᴮ) > eclipse_threshold) # in sun -> Make consistent across everything!

        if (ecl < 0.9)
            println("ERROR!: ecl: $ecl but s = ", norm(𝐬̂ᴮ), " at time $t")
            # @infiltrate
        end
        flags.in_sun = true 
        if !flags.diodes_calibrated
            ### CALIBRATE DIODES
            q, R = triad(ŝᴵ, B̂ᴵ, 𝐬̂ᴮ, Bᴮ) # Normalizes all input Vectors
            mode = diode_cal 
            cont = TRIVIAL(0.0)

            println("Making new diode: ")
            est = new_diode_calib(albedo, sens, system, q, satellite_estimate) # Init DIODE_CALIB
            est.inertial_vecs = [ŝᴵ B̂ᴵ]'
            est.body_vecs = [𝐬̂ᴮ Bᴮ]'
            est.current_meas = sens.diodes

            if norm(𝐬̂ᴮ) < 0.01
                println("Somehow still not actually in sun...")
                @infiltrate
            end

            flags.calibrating = true
            global temp_count = 0 # Ensures that everytime you enter into sun you give yourself a second to calm before finishing
        else 
            mode = mekf 
            q, R = triad(ŝᴵ, B̂ᴵ, 𝐬̂ᴮ, Bᴮ)
            if isa(current_data, DIODE_CALIB) # Finished diode cal during sunlit cycle 
                # est = new_mekf(current_data)  
                est = new_mekf(albedo, sens, system, q, satellite_estimate)
            elseif isa(current_data, MEKF)
                est = current_data
            else    # Starting out of an eclipse 
                est = new_mekf(albedo, sens, system, q, satellite_estimate)
            end

            # est = new_diode_calib(albedo, sens, system, q, satellite_estimate)
            est.inertial_vecs = [ŝᴵ B̂ᴵ]'
            est.body_vecs = [𝐬̂ᴮ Bᴮ]'
            est.current_meas = sens.diodes
            est.pos = sens.gps
            est.ang_vel = sens.gyro


            cont = TRIVIAL(0.0)
            flags.calibrating = false
        end

    # elseif ((flags.diodes_calibrated) && (norm(ω) > deg2rad(10.0))) || ((!flags.diodes_calibrated) && (norm(ω) > deg2rad(20.0)))
    #     mode = detumble
    #     cont = DETUMBLER(sens.gyro, Bᴮ, system._dt)  
    #     est = TRIVIAL(0.0)
    #     flags.detumbling = true
    else# Need to ensure diode calibration and MEKF (β) are used as initial conditions for the next loop
        flags.in_sun = false 
        mode = chill
        cont = TRIVIAL(0.0)
        est = TRIVIAL(0.0)
        flags.calibrating = false
    end

    return mode, cont, est, flags, satellite_estimate
end



# HELPER  Functions
function estimate_sun_vector(sens::SENSORS, sat_est::SATELLITE)
    if norm(sens.diodes) > eclipse_threshold  # If not eclipsed
        x₁ = (sens.diodes[1]/sat_est.diodes.calib_values[1])
        x₂ = (sens.diodes[2]/sat_est.diodes.calib_values[2])
        y₁ = (sens.diodes[3]/sat_est.diodes.calib_values[3])
        y₂ = (sens.diodes[4]/sat_est.diodes.calib_values[4])
        z₁ = (sens.diodes[5]/sat_est.diodes.calib_values[5])
        z₂ = (sens.diodes[6]/sat_est.diodes.calib_values[6])

        # sun_vec_est = [x₁ - x₂;
        #                 y₁ - y₂;
        #                 z₁ - z₂]
        # sun_vec_est = [ (sens.diodes[1]/sat_est.diodes.calib_values[1]) - (sens.diodes[2]/sat_est.diodes.calib_values[2]);
        #                 (sens.diodes[3]/sat_est.diodes.calib_values[3]) - (sens.diodes[4]/sat_est.diodes.calib_values[4]);
        #                 (sens.diodes[5]/sat_est.diodes.calib_values[5]) - (sens.diodes[6]/sat_est.diodes.calib_values[6])]


        sun_vec_est = [(x₁*cos(-pi/4) + z₁*cos(pi/4) + x₂*cos(3*pi/4) + z₂*cos(-3*pi/4));    #(z₁*cos(pi/4) + x₂*cos(3*pi/4) + z₂*cos(5*pi/4) + x₁*cos(7*pi/4))
                        y₁ - y₂;
                       (x₁*cos(3*pi/4) + z₁*cos(pi/4) + x₂*cos(-pi/4) + z₂*cos(-3*pi/4))] # -(x₁*cos(pi/4) + z₁*cos(3*pi/4) + x₂*cos(5*pi/4) + z₂*cos(7*pi/4))]    

        sun_vec_est /= norm(sun_vec_est)
    else
        # Assume that we are in eclipse
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
        (Summary)

        Arguments:
            - 

        Returns:
            -
    """
    # Method for estimating the rotation matrix between two reference frames 
    #   Relies only on a single pair of vectors in each frame
    # Inputs: 
    #   - r₁ᴵ, r₂ᴵ: Pair of vectors in the Newtonian (inertial) frame     | [3,]
    #   - r₁ᴮ, r₂ᴮ: Corresponding pair of vectors in body frame           | [3,]    
    # Outputs:
    #   - R: A directed cosine matrix (DCM) representing the rotation     | [3 x 3]
    #           between the two frames 
    #   - q: A quaternion (scalar last) representing the rotation         | [4,]
    #           between the two frames

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


# For detumbler (combine with photodiode method...?)
function check_if_finished(sens::SENSORS, thresh)
    return (norm(sens.gyro) < deg2rad(thresh))
end

# For magnetometer 
function check_if_finished(sat_est::SATELLITE)
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

# For photodiodes ( combine with detumbler)
function check_if_finished(data::DIODE_CALIB, sat_est::SATELLITE)
    return (norm(sat_est.covariance) < 0.0065) #7)
end