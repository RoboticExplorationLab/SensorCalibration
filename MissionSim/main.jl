# Main file for mission simulator

using Plots
using LinearAlgebra         
using SatelliteDynamics     # Allows for coordinate transforms, eclipse, sun position stuff, etc...
using Random, Distributions
using ProgressMeter
using EarthAlbedo
using MAT

using Infiltrator

using Random
# Random.seed!(12973)


include("mag_field.jl")            # Contains IGRF13 stuff
include("rotationFunctions.jl")    # Contains general functions for working with quaternions

include("CustomStructs.jl");         using .CustomStructs
include("system_config_file.jl"); 

include("Estimator/Estimator.jl");   using .Estimator
include("Simulator/Simulator.jl");   using .Simulator
include("Controller/Controller.jl"); using .Controller

@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, full_control, finished, chill)
Base.to_index(s::Operation_mode) = Int(s)


# Should I create a run() function that does initialization and then calls a loop/main() function  -> YES
#   that takes in arguments? Then I can just adjust the arguments...?  -> MAYBE
function main()

    println("Initializing Values")

    # Generate our satellite 
    magnetometer = MAGNETOMETER([_a, _b, _c], [_ρ, _λ, _ϕ], [_βx₀, _βy₀, _βz₀]) 
    diodes = DIODES(_sensor_scale_factors, _azi_angles, _elev_angles)
    satellite_truth = SATELLITE(_J, magnetometer, diodes) 

    operation_mode = mag_cal    # First mode: Calibrate the magnetometer

    # Initialize histories to max lenth 
    x_hist = zeros(length(x0), _max_sim_length + 1)  
    mekf_hist = zeros(7, _max_sim_length)

    # Initialize satellite estimates 
    m = MAGNETOMETER(ones(3), zeros(3), zeros(3)) # Assume no bias, no non-orthogonality, and a unit scale factor 
    m = initialize(m)
    
    
    d = DIODES(ones(_num_diodes), zeros(_num_diodes), zeros(_num_diodes))
    d = initialize(satellite_truth, d) 
    _̂J = _J
    satellite_estimate = SATELLITE(_̂J, m, d) # Jhat = J for now

    estimates_hist = Array{SATELLITE,1}(undef, _max_sim_length)
    sensors_hist = Array{SENSORS, 1}(undef, _max_sim_length)
    ground_truth_hist = Array{GROUND_TRUTH, 1}(undef, _max_sim_length)

    sun_vec_est = zeros(3, _max_sim_length)

    for i = 1:_max_sim_length 
        new_sat = SATELLITE(_̂J, 
                    MAGNETOMETER(zeros(3), zeros(3), zeros(3)),
                    DIODES(zeros(_num_diodes), zeros(_num_diodes), zeros(_num_diodes)))
        estimates_hist[i] = new_sat

        new_sens = SENSORS(zeros(3), zeros(_num_diodes), zeros(3), zeros(3))
        sensors_hist[i] = new_sens 

        new_truth = GROUND_TRUTH(_epc, zeros(Float32, 3), zeros(Float32, 3), zeros(Float32, 3))
        ground_truth_hist[i] = new_truth
    end
    

    # Generate ALBEDO data 
    #   refl_dict = matread("../../Earth_Albedo_Model/Processed_Data/tomsdata2005/2005/ga050101-051231.mat")
    refl_dict = matread("refl.mat")     # Same as ^ but saved in local directory
    refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])

    cell_centers_ecef = get_albedo_cell_centers()
    albedo = ALBEDO(refl, cell_centers_ecef)

    # Set satellite initial estimates? 
    sim = SIM(1.0)

    x_hist[:, 1] = x0
    updated_data = MAG_CALIB(0.0, 0.0)

    # operation_mode = diode_cal
    # updated_data = initialize(albedo, x0, SYSTEM)
    # satellite_estimate.magnetometer = satellite_truth.magnetometer  ###############
    # operation_mode = detumble

    temp = zeros(3, _max_sim_length) # Used for un-corrected mag field
    temp_cor = zeros(3, _max_sim_length)
    count = 0

    temp_flag = false
    progress_bar = Progress(_max_sim_length)









    # Design alternative with if/else structure instead of uniformity to compare complexity 
    for i = 1:_max_sim_length
        state = x_hist[:, i]
        t = ground_truth_hist[i].t_hist + (i - 1) * _dt

        # STATE MACHINE
        # mutable stuct FLAGS (sun_flag, mag, diode, triad?)
        # controller, estimator, operation_mode = update_operation_mode(state, flags)
        #   + TRIAD flag?   |   estimator == updated_data?   |
        # Detumble -> Calib Mag --(SUN)-> Calib Diodes --(SUN)-> MEKF 
        # if (diodes_calib) subtract off bias in detumble?

        # Add in SIM vs HARDWARE options
        truth, sensors = generate_measurements(sim, satellite_truth, albedo, state, t, CONSTANTS, _dt)

        temp[:, i] = sensors.magnetometer
        sun_vec_est[:, i] = estimate_sun_vector(sensors, satellite_estimate)
        # sun_vec_est[:, i] = truth.ŝᴮ_hist

        # if flags.mag_calibrated
        sensors.magnetometer = correct_mag_field(satellite_estimate, sensors.magnetometer) # Correct using our calibration estimate
        # DO I set sun vector = 0 and then generate here, or generate whenever i need it? Second seems more accurate...
        temp_cor[:, i] = sensors.magnetometer

        estimators = sensors_to_estimators(sensors, truth, updated_data, satellite_estimate)

        old_estimate = deepcopy(satellite_estimate)

        if operation_mode == mag_cal
            if i > Int(round(2 * orbit_period(oe0[1])))
                satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode], true)
            elseif i % 300 == 0 
                satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode], false)
            else
                estimator_finished = false
            end
        else
            satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode])
        end


        # if i % 60 == 0 
        #     satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode])
        # else 
        #     estimator_finished = false
        # end
        # estimator_finished = false

        # if i >  Int(round(2 * orbit_period(oe0[1])))# 5750
        #     estimator_finished = true 
        # end

        change = 0

        if operation_mode == diode_cal 
            change, estimator_finished = check_if_finished(old_estimate, satellite_estimate, updated_data)
            if count < 2500
                estimator_finished = false 
            end
        end
        # estimator_finished = false 

        controllers = sensors_to_controllers(sensors, updated_data, _dt)

        control_input, controller_finished = generate_command(controllers[operation_mode], state)
        if (operation_mode == detumble)
            if count < 100
                controller_finished = false
            else
                change, controller_finished = check_if_finished(x_hist, i) # UPDATE to be mag(w) < thresh || run long enough? 
            end
        end

        new_state = rk4(satellite_truth, state, control_input, t, _dt)

        # Update histories
        x_hist[:, i + 1] = new_state
        x_hist[7:10, i + 1] = x_hist[7:10, i + 1] / norm(x_hist[7:10, i + 1]) # Normalize quaternions


        if operation_mode == diode_cal 
            mekf_hist[1:4, i] = updated_data.sat_state[1:4]   # Quats
            mekf_hist[5:7, i] = updated_data.sat_state[5:7]   # Bias
        else 
            mekf_hist[:, i] = zeros(7)
        end

        sensors_hist[i] = SENSORS(sensors.magnetometer, sensors.diodes, sensors.gyro, sensors.gps)
        ground_truth_hist[i] = GROUND_TRUTH(t, truth.Bᴵ_hist, truth.sᴵ_hist, truth.ŝᴮ_hist)
        estimates_hist[i] = deepcopy(satellite_estimate)

        # Update mode - Currently ugly b/c cant iterate an @enum value... yet
        if (operation_mode == mag_cal) && (estimator_finished)
            println("\nFinished Calibrating Magnetometer")
            operation_mode = detumble # detumble
            updated_data = TRIVIAL(0.0)
            count = 0 
            report_on_magnetometer(satellite_truth, satellite_estimate, estimates_hist[1]) 

        elseif (operation_mode == detumble) && (controller_finished)
            println("\nFinished Detumbling")
            operation_mode = diode_cal
            count = 0
            # updated_data = initialize(albedo, new_state, SYSTEM)
            updated_data = new_diode_calib(albedo, sensors, SYSTEM)
            report_on_detumbler(x_hist[:, 1:i], sensors_hist)


        elseif (operation_mode == diode_cal) && (estimator_finished)
            println("\nFinished Calibrating Diodes")
            operation_mode = finished 
            report_on_diodes(satellite_truth, estimates_hist[1:i])
            # temp_flag = true

        elseif (operation_mode == full_control) && (controller_finished)
            operation_mode = finished
        end

        if operation_mode == finished
            x_hist = x_hist[:, 1:(i+1)]   # Save through most recent 
            # save_global_variables(satellite_truth, estimates_hist[1])  # Used to test magnetometer downsampling
            break
        end

        count += 1

        if norm(truth.ŝᴮ_hist) < 0.01
            ecl = true
        else
            ecl = false
        end
        sun_error = norm(truth.ŝᴮ_hist - sun_vec_est[:, i])
        notes = string("Change: $change    \t|   Ecl: $ecl   \t|   Sun vector error: $sun_error")
        ProgressMeter.next!(progress_bar; showvalues = [(:Mode, operation_mode), (:Iteration, i), (:Notes, notes)])
    end





    # temp = 1 # in case I need to return something extra for debugging
    return sun_vec_est, mekf_hist, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist
end


# FIX -> MAG_CALIB
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::MAG_CALIB, sat_est::SATELLITE)
    # Assumes that IGRF IS PERFECT
    """ [MAG_CALIB; TRIVIAL; DIODE_CALIB; MEKF] """
    # MAG_CALIB: B_meas, B_pred, current 
    B_meas = sens.magnetometer
    B_pred = truth.Bᴵ_hist
    currents = sens.diodes 
    mag_calib = MAG_CALIB(B_meas, B_pred) 

    # TRIVIAL: Garbage 
    triv = TRIVIAL(0.0)
    mekf = MEKF(1.0)

    # return [mag_calib; triv; diode_calib; mekf]
    return [mag_calib; mag_calib; mag_calib; mag_calib]
end


####################
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::DIODE_CALIB, sat_est::SATELLITE)

    # Assumes that IGRF IS PERFECT
    """ [MAG_CALIB; TRIVIAL; DIODE_CALIB; MEKF] """
    # MAG_CALIB: B_meas, B_pred, current 
    B_meas = sens.magnetometer  # IN BODY FRAME (== Bᴮ)
    B_pred = truth.Bᴵ_hist      # IN INTERTIAL FRAME! (== Bᴵ)
    mag_calib = MAG_CALIB(B_meas, B_pred) 

    sᴵ = truth.sᴵ_hist
    ŝᴮ = estimate_sun_vector(sens, sat_est)
    # ŝᴮ = truth.ŝᴮ_hist 


    # TRIVIAL: Garbage 
    triv = TRIVIAL(0.0)

    diode_calib = DIODE_CALIB(data.albedo, 
                                data.sat_state, 
                                data.covariance, 
                                [truth.sᴵ_hist truth.Bᴵ_hist]',
                                sens.gyro, 
                                [ŝᴮ sens.magnetometer]',
                                sens.diodes, 
                                data.W, 
                                data.V,
                                SYSTEM._dt, 
                                data.time,
                                length(sens.diodes),
                                sens.gps, 
                                data.first_pass)

    # MEKF: Lots, but less lots  -> Update
    mekf = MEKF(1.0)
    
    return [triv; diode_calib; diode_calib; triv]
end
# sensors_to_estimators for TRIVIAL needs to be written
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::TRIVIAL, sat_est::SATELLITE)
    triv = TRIVIAL(0.0)
    return [triv; triv; triv; triv]
end



function check_if_finished(old_est::SATELLITE, new_est::SATELLITE, data::DIODE_CALIB)

    if data.first_pass == true
        return 0.0, false
    end
    if sum(abs.(data.inertial_vecs[1,:])) > 0.01 # Cant finish in eclipse
        old_calib = old_est.diodes.calib_values 
        new_calib = new_est.diodes.calib_values 
        δcalib = norm(new_calib - old_calib)

        old_azi = old_est.diodes.azi_angles 
        new_azi = new_est.diodes.azi_angles
        δazi = norm(new_azi - old_azi)

        old_elev = old_est.diodes.elev_angles
        new_elev = new_est.diodes.elev_angles
        δelev = norm(new_elev - old_elev)

        change = δcalib + δazi + δelev

        if (change < 1e-5) #(δcalib < 1e-6) && (δazi < 1e-6) && (δelev < 1e-6)
            return change, true 
        else
            return change, false 
        end
    else
        return 0.0, false
    end
end
function check_if_finished(old_est::SATELLITE, new_est::SATELLITE, data::TRIVIAL)
    return 0.0, true
end
function check_if_finished(hist,i)
    if i < 101
        return false
    else
        # println(norm(hist[11:13, i]))
        ωs = hist[11:13, (i-100):i]
        T = size(ωs, 2)
        norms = zeros(T)
        for j = 1:T
            norms[j] = norm(ωs[:, j])
        end

        μω = sum(norms) / T 

        change = sum(abs.(norms .- μω))

        current_ang_vel = norm(hist[11:13, i])

        if current_ang_vel < deg2rad(6.0) #  change < 0.00075
            # println("Finished controller: ", sum(abs.(norms .- μω)))
            return change, true 
        else
            return change, false
        end

    end
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

function estimate_sun_vector(sens::SENSORS, sat_est::SATELLITE)
    if norm(sens.diodes) > 0.01  # If not eclipsed
        # surface_normals, measured_current = [], []
        # empty_flag = true
        # ϵ, α, C = sat_est.diodes.elev_angles, sat_est.diodes.azi_angles, sat_est.diodes.calib_values
        # for i = 1:length(sens.diodes) 
        #     if sens.diodes[i] != 0.0 
        #         if empty_flag
        #             surface_normals = [(cos(ϵ[i])*cos(α[i])) (cos(ϵ[i])*sin(α[i])) sin(ϵ[i])]
        #             measured_current = sens.diodes[i] / C[i]
        #             empty_flag = false
        #         else
        #             surface_normal = [(cos(ϵ[i])*cos(α[i])) (cos(ϵ[i])*sin(α[i])) sin(ϵ[i])]
        #             surface_normals = [surface_normals; surface_normal]
        
        #             measured_current = [measured_current; sens.diodes[i]/C[i]]
        #         end
        #     end 
        # end

        # if size(measured_current, 1) < 3
        #     println("ERROR! Not enough diodes illuminated")
        #     sun_vec_est = [0; 0; 0]
        # else
        #     sun_vec_est = ((surface_normals' * surface_normals)^(-1)) * (surface_normals' * measured_current)
        # end
        sun_vec_est = [ (sens.diodes[1]/sat_est.diodes.calib_values[1]) - (sens.diodes[2]/sat_est.diodes.calib_values[2]);
                        (sens.diodes[3]/sat_est.diodes.calib_values[3]) - (sens.diodes[4]/sat_est.diodes.calib_values[4]);
                        (sens.diodes[5]/sat_est.diodes.calib_values[5]) - (sens.diodes[6]/sat_est.diodes.calib_values[6])]

        sun_vec_est /= norm(sun_vec_est)

    else
        # Assume that we are in eclipse
        sun_vec_est = [0; 0; 0]
    end
        

    return sun_vec_est # Unit - ŝᴮ
end




# NEEDS TO BE UPDATED EVERY TIME WE ADD SOMETHING or ADJUST ORDER (maybe make a dict at top with name -> index and use for all?)
function sensors_to_controllers(sens::SENSORS, dt)
    """[Triv; detumbler; triv; mag]"""

    det = DETUMBLER(sens.gyro, sens.magnetometer, dt) 
    t = TRIVIAL(1.0)
    return [t; det; t; t]
end
# Duplicate with subtractability
function sensors_to_controllers(sens::SENSORS, data::DIODE_CALIB, dt)
    """[Triv; detumbler; triv; mag]"""

    gyro = sens.gyro - data.sat_state[5:7]
    det = DETUMBLER(gyro, sens.magnetometer, dt) 
    t = TRIVIAL(1.0)
    return [t; det; t; t]
end
function sensors_to_controllers(sens::SENSORS, data, dt)
    """[Triv; detumbler; triv; mag]"""

    gyro = sens.gyro #- data.sat_state[5:7]
    det = DETUMBLER(gyro, sens.magnetometer, dt) 
    t = TRIVIAL(1.0)
    return [t; det; t; t]
end








function report_on_magnetometer(truth::SATELLITE, est::SATELLITE, init::SATELLITE)
    
    aᶠ, bᶠ, cᶠ = est.magnetometer.scale_factors
        aᶠ, bᶠ, cᶠ = round(aᶠ, sigdigits = 3), round(bᶠ, sigdigits = 3), round(cᶠ, sigdigits = 3)
    ρᶠ, λᶠ, ϕᶠ = est.magnetometer.non_ortho_angles
        ρᶠ, λᶠ, ϕᶠ = round(rad2deg(ρᶠ), sigdigits = 3), round(rad2deg(λᶠ), sigdigits = 3), round(rad2deg(ϕᶠ), sigdigits = 3)
    βxᶠ, βyᶠ, βzᶠ = est.magnetometer.bias   
        βxᶠ, βyᶠ, βzᶠ = round(βxᶠ, sigdigits = 3), round(βyᶠ, sigdigits = 3), round(βzᶠ, sigdigits = 3)

    a, b, c = truth.magnetometer.scale_factors
        a, b, c = round(a, sigdigits = 3), round(b, sigdigits = 3), round(c, sigdigits = 3)
    ρ, λ, ϕ = truth.magnetometer.non_ortho_angles
        ρ, λ, ϕ = round(rad2deg(ρ), sigdigits = 3), round(rad2deg(λ), sigdigits = 3), round(rad2deg(ϕ), sigdigits = 3)
    βx, βy, βz = truth.magnetometer.bias
        βx, βy, βz = round(βx, sigdigits = 3), round(βy, sigdigits = 3), round(βz, sigdigits = 3)

    a₀, b₀, c₀ = init.magnetometer.scale_factors
        a₀, b₀, c₀ = round(a₀, sigdigits = 3), round(b₀, sigdigits = 3), round(c₀ , sigdigits = 3)
    ρ₀, λ₀, ϕ₀ = init.magnetometer.non_ortho_angles
        ρ₀, λ₀, ϕ₀ = round(rad2deg(ρ₀), sigdigits = 3), round(rad2deg(λ₀), sigdigits = 3), round(rad2deg(ϕ₀) , sigdigits = 3)
    βx₀, βy₀, βz₀ = init.magnetometer.bias
        βx₀, βy₀, βz₀ = round(βx₀, sigdigits = 3), round(βy₀, sigdigits = 3), round(βz₀ , sigdigits = 3)

    println("__________________________________________________________________________")
    println("___PARAM___|___Truth____|__Final Guess__|___Init Guess___|__Improved?__")
    println("     a     |   $a \t|    $aᶠ    \t|    $a₀         | ", abs(a - aᶠ) < abs(a - a₀) ? "better!" : "worse!")
    println("     b     |   $b \t|    $bᶠ    \t|    $b₀         | ", abs(b - bᶠ) < abs(b - b₀) ? "better!" : "worse!")
    println("     c     |   $c \t|    $cᶠ    \t|    $c₀         | ", abs(c - cᶠ) < abs(c - c₀) ? "better!" : "worse!")
    println("     ρ°    |   $ρ \t|    $ρᶠ    \t|    $ρ₀         | ", abs(ρ - ρᶠ) < abs(ρ - ρ₀) ? "better!" : "worse!")
    println("     λ°    |   $λ \t|    $λᶠ    \t|    $λ₀         | ", abs(λ - λᶠ) < abs(λ - λ₀) ? "better!" : "worse!")
    println("     ϕ°    |   $ϕ \t|    $ϕᶠ    \t|    $ϕ₀         | ", abs(ϕ - ϕᶠ) < abs(ϕ - ϕ₀) ? "better!" : "worse!")
    println("     βx    |   $βx \t|    $βxᶠ   \t|    $βx₀         | ", abs(βx - βxᶠ) < abs(βx - βx₀) ? "better!" : "worse!")
    println("     βy    |   $βy \t|    $βyᶠ   \t|    $βy₀         | ", abs(βy - βyᶠ) < abs(βy - βy₀) ? "better!" : "worse!")
    println("     βz    |   $βz \t|    $βzᶠ   \t|    $βz₀         | ", abs(βz - βzᶠ) < abs(βz - βz₀) ? "better!" : "worse!")
    println("__________________________________________________________________________")
end
function report_on_detumbler(x, sensors_hist::Array{SENSORS,1})

    N = size(x, 2) - 1
    ŵ = zeros(3, N)
    mag_w = zeros(N)

    for i = 1:N 
        ŵ[:, i] = rad2deg.(sensors_hist[i].gyro)
        mag_w[i] = rad2deg(norm(x[11:13, i]))
    end

    w = plot(rad2deg.(x[11:13, :])', color = [:red :blue :green], label = ["wₓ" "wy" "wz"], title = "Ang Vel (deg/sec)")
    w = plot!(ŵ', linestyle = :dash, color = [:red :blue :green], label = ["ŵₓ" "ŵy" "ŵz"])
    w = plot!(mag_w, linestyle = :dot, color = :black, label = "|w_true|")
    display(w)

    bias = x[14:16, :]

    display(plot(bias', title = "Bias (random walk)"))
    
end
function report_on_diodes(truth::SATELLITE, est_hist::Array{SATELLITE,1})
    N = size(est_hist, 1)
    c_est, α_est, ϵ_est = zeros(6, N), zeros(6, N), zeros(6, N)

    for i = 1:N
        c_est[:, i] = est_hist[i].diodes.calib_values
        α_est[:, i] = rad2deg.(est_hist[i].diodes.azi_angles)
        ϵ_est[:, i] = rad2deg.(est_hist[i].diodes.elev_angles)
    end
    ϵ₀ = rad2deg.([0.0; 0.0;  0.0;    0.0;    (pi/4); (-pi/4)])
    α₀ = rad2deg.([0.0; 0.0;  0.0;    0.0;    (pi/4); (-pi/4)])
    c₀ = ones(6)


    # c₀, ϵ₀, α₀ = c_est[:, 1], ϵ_est[:, 1], α_est[:, 1]
    cᶠ, ϵᶠ, αᶠ = c_est[:, end], ϵ_est[:, end], α_est[:, end]
    c, ϵ, α = truth.diodes.calib_values, rad2deg.(truth.diodes.elev_angles), rad2deg.(truth.diodes.azi_angles)

    println("\n----- DIODE REPORT -----")
    for i = 1:_num_diodes
        print("Diode $i:")
        print( abs(c₀[i] - c[i]) ≤ abs(cᶠ[i] - c[i]) ? "C: Worse!" : "C: Better!"); print("  (", round(c[i] ,sigdigits = 3), " & ", round(cᶠ[i], sigdigits = 3), " vs ", round(c₀[i], sigdigits = 3) ," (init))    \t| ")
        print( abs(ϵ₀[i] - ϵ[i]) ≤ abs(ϵᶠ[i] - ϵ[i]) ? "E: Worse! " : "E: Better!"); print("  (", round(ϵ[i] ,sigdigits = 3), " & ", round(ϵᶠ[i], sigdigits = 3), " vs ", round(ϵ₀[i], sigdigits = 3) ," (init)) \t| ")
        print( abs(α₀[i] - α[i]) ≤ abs(αᶠ[i] - α[i]) ? "A: Worse! " : "A: Better!"); print("  (", round(α[i] ,sigdigits = 3), " & ", round(αᶠ[i], sigdigits = 3), " vs ", round(α₀[i], sigdigits = 3) ," (init))\n")
    end
    println("------------------------")

    c_off = 0.2
    e_off = 3.0
    a_off = 3.0
    cp, ap, ep = Array{Plots.Plot{Plots.GRBackend}, 1}(undef, _num_diodes), Array{Plots.Plot{Plots.GRBackend}, 1}(undef, _num_diodes), Array{Plots.Plot{Plots.GRBackend}, 1}(undef, _num_diodes)
    for i = 1:_num_diodes
        cp[i] = plot(c_est[i, :], title = "Scale Factor (C)", label = false); cp[i] = hline!([c₀[i]], linestyle = :dot, label = false); cp[i] = hline!([c[i]], ylim = [c[i] - c_off, c[i] + c_off], linestyle = :dash, label = false)
        ep[i] = plot(ϵ_est[i, :], title = "Elevation Angle (ϵ)", label = false); ep[i] = hline!([ϵ₀[i]], linestyle = :dot, label = false); ep[i] = hline!([ϵ[i]], ylim = [ϵ[i] - e_off, ϵ[i] + e_off], linestyle = :dash, label = false); 
        ap[i] = plot(α_est[i, :], title = "Azimuth Angle (α)", label = false); ap[i] = hline!([α₀[i]], linestyle = :dot, label = false); ap[i] = hline!([α[i]], ylim = [α[i] - a_off, α[i] + a_off], linestyle = :dash, label = false); 
    end
    cp[2].series_list[1][:label] = "Estimate"; cp[2].series_list[2][:label] = "Initial Guess"; cp[2].series_list[3][:label] = "Truth"
    ep[2].series_list[1][:label] = "Estimate"; ep[2].series_list[2][:label] = "Initial Guess"; ep[2].series_list[3][:label] = "Truth"
    ap[2].series_list[1][:label] = "Estimate"; ap[2].series_list[2][:label] = "Initial Guess"; ap[2].series_list[3][:label] = "Truth"

    display(plot(cp..., layout = (3, 2)))
    display(plot(ep..., layout = (3, 2)))
    display(plot(ap..., layout = (3, 2)))
end
function report_on_mekf(mekf_hist, state) 
    N = size(mekf_hist, 1)
    poses, biases = zeros(4, N), zeros(3, N)

    for i = 1:N
        poses[:, i] = state[7:10, i]
        biases[:,i] = state[14:16,i]
    end
end



sun_vec_est, mekf_hist, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist = main();

display(plot(x_hist[1:3, :]', title = "Position"))


N = size(x_hist, 2) - 1
Bᴵ = zeros(3, N)
sᴵ = zeros(3, N)
for i = 1:N
    Bᴵ[:, i] = ground_truth_hist[i].Bᴵ_hist
    sᴵ[:, i] = ground_truth_hist[i].sᴵ_hist
end

tru1 = plot(Bᴵ', title = "Bᴵ")
tru2 = plot(sᴵ', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1)))

# report_on_diodes(satellite_truth, estimates_hist)

##### MAGNETOMETER CALIBRATION
N = size(x_hist, 2) - 1
Bᴵ = zeros(3, N)
sᴵ = zeros(3, N)

diodes = zeros(6, N)
gyro = zeros(3, N)
Bᴮ = zeros(3, N)
ŝᴮ = zeros(3, N)
currents = zeros(6, N)
for i = 1:N
    Bᴮ[:, i] = sensors_hist[i].magnetometer
    gyro[:, i] = sensors_hist[i].gyro 
    diodes[:, i] = sensors_hist[i].diodes
    Bᴵ[:, i] = ground_truth_hist[i].Bᴵ_hist
    sᴵ[:, i] = ground_truth_hist[i].sᴵ_hist
    ŝᴮ[:, i] = ground_truth_hist[i].ŝᴮ_hist
    currents[:, i] = sensors_hist[i].diodes
end
tru1 = plot(Bᴵ', title = "Bᴵ")
tru2 = plot(sᴵ', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1)))

sen1 = plot(Bᴮ', title = "Bᴮ")
sen2 = plot(ŝᴮ', title = "sᴮ")
sen3 = plot(diodes', title = "Diodes")
sen4 = plot(gyro', title = "Gyro")
display(plot(sen1, sen2, sen3, sen4, layout = (4,1)))


a, b, c = satellite_truth.magnetometer.scale_factors 
ρ, λ, ϕ = satellite_truth.magnetometer.non_ortho_angles
a_est, b_est, c_est = zeros(N), zeros(N), zeros(N)
ρ_est, λ_est, ϕ_est = zeros(N), zeros(N), zeros(N)
for i = 1:N
    a_est[i], b_est[i], c_est[i] = estimates_hist[i].magnetometer.scale_factors
    ρ_est[i], λ_est[i], ϕ_est[i] = estimates_hist[i].magnetometer.non_ortho_angles
end

# ap = plot(a_est, title = "ã", label = "Estimate"); ap = hline!([a_est[1]], linestyle = :dot, label = "Initial Guess"); ap = hline!([a], linestyle = :dash, label = "Truth", ylims = [0.5, 1.5])
# bp = plot(b_est, title = "b̃", label = "Estimate"); bp = hline!([b_est[1]], linestyle = :dot, label = "Initial Guess"); bp = hline!([b], linestyle = :dash, label = "Truth", ylims = [0.5, 1.5])
# cp = plot(c_est, title = "c̃", label = "Estimate"); cp = hline!([c_est[1]], linestyle = :dot, label = "Initial Guess"); cp = hline!([c], linestyle = :dash, label = "Truth", ylims = [0.5, 1.5])
# display(plot(ap, bp, cp, layout = (3, 1)))


# ρp = plot(rad2deg.(ρ_est), title = "ρ", label = "Estimate"); ρp = hline!([rad2deg(ρ_est[1])], linestyle = :dot, label = "Initial Guess"); ρp = hline!([rad2deg(ρ)], ylims = [rad2deg(ρ)-10, rad2deg(ρ)+10], linestyle = :dash, label = "Truth")
# λp = plot(rad2deg.(λ_est), title = "λ", label = "Estimate"); λp = hline!([rad2deg(λ_est[1])], linestyle = :dot, label = "Initial Guess"); λp = hline!([rad2deg(λ)], ylims = [rad2deg(λ)-10, rad2deg(λ)+10], linestyle = :dash, label = "Truth")
# ϕp = plot(rad2deg.(ϕ_est), title = "ϕ", label = "Estimate"); ϕp = hline!([rad2deg(ϕ_est[1])], linestyle = :dot, label = "Initial Guess"); ϕp = hline!([rad2deg(ϕ)], ylims = [rad2deg(ϕ)-10, rad2deg(ϕ)+10], linestyle = :dash, label = "Truth")
# display(plot(ρp, λp, ϕp, layout = (3, 1)))



#
if N == _max_sim_length # Diodes didn't finish
    report_on_diodes(satellite_truth, estimates_hist)
end

sp = plot(sun_vec_est[:,1:N]', title = "Sun Vector (body) Estimate", label = ["x̂" "ŷ" "ẑ"])
sp = plot!(ŝᴮ[:, 1:N]', linestyle = :dash, label = ["x" "y" "z"])
display(plot(sp))

err = zeros(3, N);
n_err = zeros(N);

for i = 1:N
    err[:, i] = ŝᴮ[:, i] - sun_vec_est[:, i];
    n_err[i] = norm(err[:, i]);
end

ep = plot(err', title = "Error in Sun vector estimate", label = ["x" "y" "z"])
ep2 = plot(n_err, title = "Norm of Error")
display(plot(ep, ep2, layout = (2,1)))

qp1 = plot(mekf_hist[1:4, 1:N]', title = "Attitude Estimate"); qp2 = plot(x_hist[7:10,:]', title = "Truth");  qp3 = plot( (x_hist[7:10,2:(end)] - mekf_hist[1:4, 1:N])', title = "Dif");
display(plot(qp1, qp2, qp3, layout = (3,1)))

bp1 = plot(mekf_hist[5:7, 1:N]',   title = "Bias Estimate"  ); bp2 = plot(x_hist[14:16,:]', title = "Truth"); bp3 = plot( (x_hist[14:16,2:(end)] - mekf_hist[5:7, 1:N])', title = "Dif", ylim = [-0.02, 0.02]);
display(plot(bp1, bp2, bp3, layout = (3,1)))


d1 = plot(currents[1,:], title = "Diode Currents")
d2 = plot(currents[2,:], title = "Diode Currents")
d3 = plot(currents[3,:], title = "Diode Currents")
d4 = plot(currents[4,:], title = "Diode Currents")
d5 = plot(currents[5,:], title = "Diode Currents")
d6 = plot(currents[6,:], title = "Diode Currents")
display(plot(d1, d2, d3, d4, d5, d6, layout = (3,2)))