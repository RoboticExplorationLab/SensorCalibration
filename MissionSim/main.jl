# Main file for mission simulator

using Plots
using LinearAlgebra         
using SatelliteDynamics     # Allows for coordinate transforms, eclipse, sun position stuff, etc...
using Random, Distributions
using ProgressMeter
using EarthAlbedo
using MAT


using Random
# Random.seed!(123414)


include("mag_field.jl")            # Contains IGRF13 stuff
include("rotationFunctions.jl")    # Contains general functions for working with quaternions

include("CustomStructs.jl");         using .CustomStructs
include("system_config_file.jl"); 
include("Estimator/Estimator.jl");   using .Estimator
include("Simulator/Simulator.jl");   using .Simulator
include("Controller/Controller.jl"); using .Controller


# Shared functions (diode_albedo, 
@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, full_control, finished)
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

    # Initialize satellite estimates 
    m = MAGNETOMETER(ones(3), zeros(3), zeros(3)) # Assume no bias, no non-orthogonality, and a unit scale factor 
    m = initialize(m)
    
    
    d = DIODES(zeros(_num_diodes), zeros(_num_diodes), zeros(_num_diodes))
    d = initialize(satellite_truth, d) 
    _̂J = _J
    satellite_estimate = SATELLITE(_̂J, m, d) # Jhat = J for now

    # satellite_initial = deepcopy(satellite_estimate)


    println("Pre-allocating history arrays")
    # Record
    estimates_hist = Array{SATELLITE,1}(undef, _max_sim_length)
    sensors_hist = Array{SENSORS, 1}(undef, _max_sim_length)
    ground_truth_hist = Array{GROUND_TRUTH, 1}(undef, _max_sim_length)
    for i = 1:_max_sim_length 
        new_sat = SATELLITE(_̂J, 
                    MAGNETOMETER(zeros(3), zeros(3), zeros(3)),
                    DIODES(zeros(_num_diodes), zeros(_num_diodes), zeros(_num_diodes)))
        estimates_hist[i] = new_sat

        new_sens = SENSORS(zeros(3), zeros(3), zeros(_num_diodes), zeros(3), zeros(3))
        sensors_hist[i] = new_sens 

        new_truth = GROUND_TRUTH(_epc, zeros(Float32, 3), zeros(Float32, 3))
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








    @showprogress "Working" for i = 1:_max_sim_length
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

        # if flags.mag_calibrated
        sensors.magnetometer = correct_mag_field(satellite_estimate, sensors.magnetometer) # Correct using our calibration estimate
        # DO I set sun vector = 0 and then generate here, or generate whenever i need it? Second seems more accurate...
        temp_cor[:, i] = sensors.magnetometer

        estimators = sensors_to_estimators(sensors, truth, updated_data)

        old_estimate = deepcopy(satellite_estimate)

        satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode])

        if operation_mode == diode_cal 
            estimator_finished = check_if_finished(old_estimate, satellite_estimate, updated_data)
            if count < 3000
                estimator_finished = false 
            end
        end
        # estimator_finished = false 

        controllers = sensors_to_controllers(sensors, _dt)

        control_input, controller_finished = generate_command(controllers[operation_mode], state)
        if (operation_mode == detumble)
            if count < 100
                controller_finished = false
            else
                controller_finished = check_if_finished(x_hist, i) # UPDATE to be mag(w) < thresh || run long enough? 
            end
        end

        new_state = rk4(satellite_truth, state, control_input, t, _dt)

        # Update histories
        x_hist[:, i + 1] = new_state
        x_hist[7:10, i + 1] = x_hist[7:10, i + 1] / norm(x_hist[7:10, i + 1]) # Normalize quaternions

        sensors_hist[i] = SENSORS(sensors.magnetometer, sensors.sun, sensors.diodes, sensors.gyro, sensors.gps)
        ground_truth_hist[i] = GROUND_TRUTH(t, truth.Bᴵ_hist, truth.sᴵ_hist)
        estimates_hist[i] = deepcopy(satellite_estimate)



        # Update mode - Currently ugly b/c cant iterate an @enum value... yet
        if (operation_mode == mag_cal) && (estimator_finished)
            println("\nFinished Calibrating Magnetometer")
            operation_mode = detumble # detumble
            updated_data = TRIVIAL(0.0)
            count = 0 
            report_on_magnetometer(satellite_truth, satellite_estimate, estimates_hist[1]) # satellite_initial)


            # satellite_estimate = satellite_truth                                ############################################
            # println("Setting estimate to truth for testing")

        elseif (operation_mode == detumble) && (controller_finished)
            println("\nFinished Detumbling")
            operation_mode = diode_cal 
            count = 0
            updated_data = initialize(albedo, new_state, SYSTEM)
            report_on_detumbler(x_hist[:, 1:i], sensors_hist)

        elseif (operation_mode == diode_cal) && (estimator_finished)
            println("\nFinished Calibrating Diodes")
            operation_mode = finished 
            report_on_diodes(satellite_truth, estimates_hist[1:i])

        elseif (operation_mode == full_control) && (controller_finished)
            operation_mode = finished
        end


        if operation_mode == finished
            x_hist = x_hist[:, 1:i]
            break
        end

        count += 1
    end

    # temp = 1 # in case I need to return something extra for debugging
    return temp, temp_cor, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist
end


# (CAN I JUST LEAVE data AMBIGUOUS AND USE AN 'if' STATEMENT? Or pass in "finished" as well and then dont return an array?)
# FIX -> MAG_CALIB
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::MAG_CALIB)
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

# Do I need to return a set if I am overloading? DIODE_CALIB
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::DIODE_CALIB)

    # Assumes that IGRF IS PERFECT
    """ [MAG_CALIB; TRIVIAL; DIODE_CALIB; MEKF] """
    # MAG_CALIB: B_meas, B_pred, current 
    B_meas = sens.magnetometer
    B_pred = truth.Bᴵ_hist
    mag_calib = MAG_CALIB(B_meas, B_pred) 

    # TRIVIAL: Garbage 
    triv = TRIVIAL(0.0)

    diode_calib = DIODE_CALIB(data.albedo, 
                                data.sat_state, 
                                data.covariance, 
                                [truth.sᴵ_hist truth.Bᴵ_hist]',
                                sens.gyro, 
                                [sens.sun sens.magnetometer]',
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
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::TRIVIAL)
    triv = TRIVIAL(0.0)
    return [triv; triv; triv; triv]
end

function check_if_finished(old_est::SATELLITE, new_est::SATELLITE, data::DIODE_CALIB)

    if data.first_pass == true
        return false
    end
    old_calib = old_est.diodes.calib_values 
    new_calib = new_est.diodes.calib_values 
    δcalib = norm(new_calib - old_calib)

    old_azi = old_est.diodes.azi_angles 
    new_azi = new_est.diodes.azi_angles
    δazi = norm(new_azi - old_azi)

    old_elev = old_est.diodes.elev_angles
    new_elev = new_est.diodes.elev_angles
    δelev = norm(new_elev - old_elev)

    # println("$δcalib $δazi $δelev")
    if (δcalib < 1e-6) && (δazi < 1e-6) && (δelev < 1e-6)
        # println("Finished Calibrating Diodes")
        return true 
    else
        return false 
    end
end

function check_if_finished(old_est::SATELLITE, new_est::SATELLITE, data::TRIVIAL)
    return true
end

function check_if_finished(hist,i)
    if i < 101
        return false
    else
        ωs = hist[11:13, (i-100):i]
        T = size(ωs, 2)
        norms = zeros(T)
        for j = 1:T
            norms[j] = norm(ωs[:, j])
        end

        μω = sum(norms) / T 

        if sum(abs.(norms .- μω)) < 0.00075
            # println("Finished controller: ", sum(abs.(norms .- μω)))
            return true 
        else
            return false
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

    return B̂ᴮ # Corrected mag field in body frame
end





# NEEDS TO BE UPDATED EVERY TIME WE ADD SOMETHING or ADJUST ORDER (maybe make a dict at top with name -> index and use for all?)
function sensors_to_controllers(sens::SENSORS, dt)
    """[Triv; detumbler; triv; mag]"""

    det = DETUMBLER(sens.gyro, sens.magnetometer, dt) 
    t = TRIVIAL(1.0)
    return [t; det; t; t]
end









function report_on_magnetometer(truth::SATELLITE, est::SATELLITE, init::SATELLITE)
    
    aᶠ, bᶠ, cᶠ = est.magnetometer.scale_factors
    ρᶠ, λᶠ, ϕᶠ = est.magnetometer.non_ortho_angles
    βxᶠ, βyᶠ, βzᶠ = est.magnetometer.bias

    a, b, c = truth.magnetometer.scale_factors
    ρ, λ, ϕ = truth.magnetometer.non_ortho_angles
    βx, βy, βz = truth.magnetometer.bias

    a₀, b₀, c₀ = init.magnetometer.scale_factors
    ρ₀, λ₀, ϕ₀ = init.magnetometer.non_ortho_angles
    βx₀, βy₀, βz₀ = init.magnetometer.bias

    println("a: |$a - $aᶠ| = ", a - aᶠ, " (vs ", a - a₀ , ") ", abs(a - aᶠ) < abs(a - a₀) ? "better!" : "worse!")
    println("b: |$b - $bᶠ| = ", b - bᶠ, " (vs ", b - b₀ , ") ", abs(b - bᶠ) < abs(b - b₀) ? "better!" : "worse!")
    println("c: |$c - $cᶠ| = ", c - cᶠ, " (vs ", c - c₀ , ") ", abs(c - cᶠ) < abs(c - c₀) ? "better!" : "worse!")
    println("ρ: |$ρ - $ρᶠ| = ", rad2deg(ρ - ρᶠ), "° (vs ", rad2deg(ρ - ρ₀) , "°) ", abs(ρ - ρᶠ) < abs(ρ - ρ₀) ? "better!" : "worse!")
    println("λ: |$λ - $λᶠ| = ", rad2deg(λ - λᶠ), "° (vs ", rad2deg(λ - λ₀) , "°) ", abs(λ - λᶠ) < abs(λ - λ₀) ? "better!" : "worse!")
    println("ϕ: |$ϕ - $ϕᶠ| = ", rad2deg(ϕ - ϕᶠ), "° (vs ", rad2deg(ϕ - ϕ₀) , "°) ", abs(ϕ - ϕᶠ) < abs(ϕ - ϕ₀) ? "better!" : "worse!")

    println("βx: |$βx - $βxᶠ| = ", βx - βxᶠ, " (vs ", βx - βx₀ , ") ", abs(βx - βxᶠ) < abs(βx - βx₀) ? "better!" : "worse!")
    println("βy: |$βy - $βyᶠ| = ", βy - βyᶠ, " (vs ", βy - βy₀ , ") ", abs(βy - βyᶠ) < abs(βy - βy₀) ? "better!" : "worse!")
    println("βz: |$βz - $βzᶠ| = ", βz - βzᶠ, " (vs ", βz - βz₀ , ") ", abs(βz - βzᶠ) < abs(βz - βz₀) ? "better!" : "worse!")
end

function report_on_detumbler(x, sensors_hist::Array{SENSORS,1})

    N = size(x, 2) - 1
    ŵ = zeros(3, N)
    mag_w = zeros(N)

    for i = 1:N 
        ŵ[:, i] = rad2deg.(sensors_hist[i].gyro)
        mag_w[i] = rad2deg(norm(x[11:13, i]))
    end

    w = plot(rad2deg.(x[11:13, :])', color = [:red :blue :green], label = ["wₓ" "wy" "wz"], title = "Ang Vel (rad/sec)")
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

    c₀, ϵ₀, α₀ = c_est[:, 1], ϵ_est[:, 1], α_est[:, 1]
    cᶠ, ϵᶠ, αᶠ = c_est[:, end], ϵ_est[:, end], α_est[:, end]
    c, ϵ, α = truth.diodes.calib_values, rad2deg.(truth.diodes.elev_angles), rad2deg.(truth.diodes.azi_angles)

    println("\n----- DIODE REPORT -----")
    println( abs(c₀[1] - c[1]) < abs(cᶠ[1] - c[1]) ? "C1: Worse!" : "C1: Better!" )
    println( abs(c₀[2] - c[2]) < abs(cᶠ[2] - c[2]) ? "C2: Worse!" : "C2: Better!" )
    println( abs(c₀[3] - c[3]) < abs(cᶠ[3] - c[3]) ? "C3: Worse!" : "C3: Better!" )
    println( abs(c₀[4] - c[4]) < abs(cᶠ[4] - c[4]) ? "C4: Worse!" : "C4: Better!" )
    println( abs(c₀[5] - c[5]) < abs(cᶠ[5] - c[5]) ? "C5: Worse!" : "C5: Better!" )
    println( abs(c₀[6] - c[6]) < abs(cᶠ[6] - c[6]) ? "C6: Worse!" : "C6: Better!" )
    println("------------------------")

    c_off = 0.25
    c₁ = plot(c_est[1, :], ylabel = "Scale Factor Estimate", xlabel = "Time", title = "Scale Factors", label = "Estimate"); c₁ = hline!([c[1]], ylim = [c[1] - c_off, c[1] + c_off], linestyle = :dash, label = "Truth"); c₁ = hline!([c₀[1]], linestyle = :dot, label = "Initial Guess")
    c₂ = plot(c_est[2, :]); c₂ = hline!([c[2]], ylim = [c[2] - c_off, c[2] + c_off], linestyle = :dash); c₂ = hline!([c₀[2]], linestyle = :dot)
    c₃ = plot(c_est[3, :]); c₃ = hline!([c[3]], ylim = [c[3] - c_off, c[3] + c_off], linestyle = :dash); c₃ = hline!([c₀[3]], linestyle = :dot)
    c₄ = plot(c_est[4, :]); c₄ = hline!([c[4]], ylim = [c[4] - c_off, c[4] + c_off], linestyle = :dash); c₄ = hline!([c₀[4]], linestyle = :dot)
    c₅ = plot(c_est[5, :]); c₅ = hline!([c[5]], ylim = [c[5] - c_off, c[5] + c_off], linestyle = :dash); c₅ = hline!([c₀[5]], linestyle = :dot)
    c₆ = plot(c_est[6, :]); c₆ = hline!([c[6]], ylim = [c[6] - c_off, c[6] + c_off], linestyle = :dash); c₆ = hline!([c₀[6]], linestyle = :dot)
    display(plot(c₁, c₂, c₃, c₄, c₅, c₆, layout = (3, 2)))


    println( abs(ϵ₀[1] - ϵ[1]) < abs(ϵᶠ[1] - ϵ[1]) ? "E1: Worse!" : "E1: Better!" )
    println( abs(ϵ₀[2] - ϵ[2]) < abs(ϵᶠ[2] - ϵ[2]) ? "E2: Worse!" : "E2: Better!" )
    println( abs(ϵ₀[3] - ϵ[3]) < abs(ϵᶠ[3] - ϵ[3]) ? "E3: Worse!" : "E3: Better!" )
    println( abs(ϵ₀[4] - ϵ[4]) < abs(ϵᶠ[4] - ϵ[4]) ? "E4: Worse!" : "E4: Better!" )
    println( abs(ϵ₀[5] - ϵ[5]) < abs(ϵᶠ[5] - ϵ[5]) ? "E5: Worse!" : "E5: Better!" )
    println( abs(ϵ₀[6] - ϵ[6]) < abs(ϵᶠ[6] - ϵ[6]) ? "E6: Worse!" : "E6: Better!" )
    println("------------------------")

    e_off = 10.0
    e₁ = plot(ϵ_est[1, :], ylabel = "Elev Estimate", xlabel = "Time", title = "Elevation Angles"); a₁ = hline!([(ϵ[1])], ylim = [(ϵ[1]) - e_off, (ϵ[1]) + e_off],  linestyle = :dash); e₁ = hline!([ϵ₀[1]], linestyle = :dot)
    e₂ = plot(ϵ_est[2, :]); e₂ = hline!([(ϵ[2])], ylim = [(ϵ[2]) - e_off, (ϵ[2]) + e_off], linestyle = :dash); e₂ = hline!([ϵ₀[2]], linestyle = :dot)
    e₃ = plot(ϵ_est[3, :]); e₃ = hline!([(ϵ[3])], ylim = [(ϵ[3]) - e_off, (ϵ[3]) + e_off], linestyle = :dash); e₃ = hline!([ϵ₀[3]], linestyle = :dot)
    e₄ = plot(ϵ_est[4, :]); e₄ = hline!([(ϵ[4])], ylim = [(ϵ[4]) - e_off, (ϵ[4]) + e_off], linestyle = :dash); e₄ = hline!([ϵ₀[4]], linestyle = :dot)
    e₅ = plot(ϵ_est[5, :]); e₅ = hline!([(ϵ[5])], ylim = [(ϵ[5]) - e_off, (ϵ[5]) + e_off], linestyle = :dash); e₅ = hline!([ϵ₀[5]], linestyle = :dot)
    e₆ = plot(ϵ_est[6, :]); e₆ = hline!([(ϵ[6])], ylim = [(ϵ[6]) - e_off, (ϵ[6]) + e_off], linestyle = :dash); e₆ = hline!([ϵ₀[6]], linestyle = :dot)
    display(plot(e₁, e₂, e₃, e₄, e₅, e₆, layout = (3, 2)))


    println( abs(α₀[1] - α[1]) < abs(αᶠ[1] - α[1]) ? "A1: Worse!" : "A1: Better!" )
    println( abs(α₀[2] - α[2]) < abs(αᶠ[2] - α[2]) ? "A2: Worse!" : "A2: Better!" )
    println( abs(α₀[3] - α[3]) < abs(αᶠ[3] - α[3]) ? "A3: Worse!" : "A3: Better!" )
    println( abs(α₀[4] - α[4]) < abs(αᶠ[4] - α[4]) ? "A4: Worse!" : "A4: Better!" )
    println( abs(α₀[5] - α[5]) < abs(αᶠ[5] - α[5]) ? "A5: Worse!" : "A5: Better!" )
    println( abs(α₀[6] - α[6]) < abs(αᶠ[6] - α[6]) ? "A6: Worse!" : "A6: Better!" )
    println("------------------------")

    a_off = 10.0
    a₁ = plot(α_est[1, :], ylabel = "Azi Estimate", xlabel = "Time", title = "Azimuth Angles"); a₁ = hline!([(α[1])], ylim = [(α[1]) - a_off, (α[1]) + a_off],  linestyle = :dash); a₁ = hline!([α₀[1]], linestyle = :dot)
    a₂ = plot(α_est[2, :]); a₂ = hline!([(α[2])], ylim = [(α[2]) - a_off, (α[2]) + a_off], linestyle = :dash); a₂ = hline!([α₀[2]], linestyle = :dot)
    a₃ = plot(α_est[3, :]); a₃ = hline!([(α[3])], ylim = [(α[3]) - a_off, (α[3]) + a_off], linestyle = :dash); a₃ = hline!([α₀[3]], linestyle = :dot)
    a₄ = plot(α_est[4, :]); a₄ = hline!([(α[4])], ylim = [(α[4]) - a_off, (α[4]) + a_off], linestyle = :dash); a₄ = hline!([α₀[4]], linestyle = :dot)
    a₅ = plot(α_est[5, :]); a₅ = hline!([(α[5])], ylim = [(α[5]) - a_off, (α[5]) + a_off], linestyle = :dash); a₅ = hline!([α₀[5]], linestyle = :dot)
    a₆ = plot(α_est[6, :]); a₆ = hline!([(α[6])], ylim = [(α[6]) - a_off, (α[6]) + a_off], linestyle = :dash); a₆ = hline!([α₀[6]], linestyle = :dot)
    display(plot(a₁, a₂, a₃, a₄, a₅, a₆, layout = (3, 2)))

end




temp, temp_cor, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist = main();

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
sᴮ = zeros(3, N)
for i = 1:N
    Bᴮ[:, i] = sensors_hist[i].magnetometer
    sᴮ[:, i] = sensors_hist[i].sun
    gyro[:, i] = sensors_hist[i].gyro 
    diodes[:, i] = sensors_hist[i].diodes
    Bᴵ[:, i] = ground_truth_hist[i].Bᴵ_hist
    sᴵ[:, i] = ground_truth_hist[i].sᴵ_hist
end
tru1 = plot(Bᴵ', title = "Bᴵ")
tru2 = plot(sᴵ', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1)))

sen1 = plot(Bᴮ', title = "Bᴮ")
sen2 = plot(sᴮ', title = "sᴮ")
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

ap = plot(a_est, title = "ã", label = "Estimate"); ap = hline!([a_est[1]], linestyle = :dot, label = "Initial Guess"); ap = hline!([a], linestyle = :dash, label = "Truth", ylims = [0.5, 1.5])
bp = plot(b_est, title = "b̃", label = "Estimate"); bp = hline!([b_est[1]], linestyle = :dot, label = "Initial Guess"); bp = hline!([b], linestyle = :dash, label = "Truth", ylims = [0.5, 1.5])
cp = plot(c_est, title = "c̃", label = "Estimate"); cp = hline!([c_est[1]], linestyle = :dot, label = "Initial Guess"); cp = hline!([c], linestyle = :dash, label = "Truth", ylims = [0.5, 1.5])
display(plot(ap, bp, cp, layout = (3, 1)))


ρp = plot(rad2deg.(ρ_est), title = "ρ", label = "Estimate"); ρp = hline!([rad2deg(ρ_est[1])], linestyle = :dot, label = "Initial Guess"); ρp = hline!([rad2deg(ρ)], ylims = [rad2deg(ρ)-10, rad2deg(ρ)+10], linestyle = :dash, label = "Truth")
λp = plot(rad2deg.(λ_est), title = "λ", label = "Estimate"); λp = hline!([rad2deg(λ_est[1])], linestyle = :dot, label = "Initial Guess"); λp = hline!([rad2deg(λ)], ylims = [rad2deg(λ)-10, rad2deg(λ)+10], linestyle = :dash, label = "Truth")
ϕp = plot(rad2deg.(ϕ_est), title = "ϕ", label = "Estimate"); ϕp = hline!([rad2deg(ϕ_est[1])], linestyle = :dot, label = "Initial Guess"); ϕp = hline!([rad2deg(ϕ)], ylims = [rad2deg(ϕ)-10, rad2deg(ϕ)+10], linestyle = :dash, label = "Truth")
display(plot(ρp, λp, ϕp, layout = (3, 1)))

# aᶠ = a_est[end]
# bᶠ = b_est[end]
# cᶠ = c_est[end]
# ρᶠ = ρ_est[end]
# λᶠ = λ_est[end]
# ϕᶠ = ϕ_est[end]
# println("a: |$a - $aᶠ| = ", a - aᶠ, " (vs ", a - a_est[1] , ") ", abs(a - aᶠ) < abs(a - a_est[1]) ? "better!" : "worse!")
# println("b: |$b - $bᶠ| = ", b - bᶠ, " (vs ", b - b_est[1] , ") ", abs(b - bᶠ) < abs(b - b_est[1]) ? "better!" : "worse!")
# println("c: |$c - $cᶠ| = ", c - cᶠ, " (vs ", c - c_est[1] , ") ", abs(c - cᶠ) < abs(c - c_est[1]) ? "better!" : "worse!")
# println("ρ: |$ρ - $ρᶠ| = ", rad2deg(ρ - ρᶠ), "° (vs ", rad2deg(ρ - ρ_est[1]) , "°) ", abs(ρ - ρᶠ) < abs(ρ - ρ_est[1]) ? "better!" : "worse!")
# println("λ: |$λ - $λᶠ| = ", rad2deg(λ - λᶠ), "° (vs ", rad2deg(λ - λ_est[1]) , "°) ", abs(λ - λᶠ) < abs(λ - λ_est[1]) ? "better!" : "worse!")
# println("ϕ: |$ϕ - $ϕᶠ| = ", rad2deg(ϕ - ϕᶠ), "° (vs ", rad2deg(ϕ - ϕ_est[1]) , "°) ", abs(ϕ - ϕᶠ) < abs(ϕ - ϕ_est[1]) ? "better!" : "worse!")

# T = [a          0               0;
#      b*sin(ρ)   b*cos(ρ)        0; 
#      c*sin(λ)   c*sin(ϕ)*cos(λ) c*cos(ϕ)*cos(λ)  ]
     

# Tᶠ = [aᶠ           0                  0;
#      bᶠ*sin(ρᶠ)   bᶠ*cos(ρᶠ)         0; 
#      cᶠ*sin(λᶠ)   cᶠ*sin(ϕᶠ)*cos(λᶠ) cᶠ*cos(ϕᶠ)*cos(λᶠ)  ]

# println(T)
# println(Tᶠ)