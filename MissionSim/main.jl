# Main file for mission simulator

# USING and INCLUDES that are used in ALL 
using Plots
using LinearAlgebra   # GENERIC/common ones
using SatelliteDynamics   # Allows for coordinate transforms, eclipse, sun position stuff, etc...
using Random, Distributions
using ProgressMeter
using EarthAlbedo
using MAT


include("mag_field.jl")            # Contains IGRF13 stuff
include("rotationFunctions.jl")    # Contains general functions for working with quaternions

include("CustomStructs.jl"); using .CustomStructs
include("system_config_file.jl"); 
include("Estimator/Estimator.jl");   using .Estimator
include("Simulator/Simulator.jl");   using .Simulator
include("Controller/Controller.jl"); using .Controller


# Shared functions (diode_albedo, 
@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, full_control, finished)
Base.to_index(s::Operation_mode) = Int(s)


function main()

    # Generate our satellite 
    magnetometer = MAGNETOMETER([_a, _b, _c], [_ρ, _λ, _ϕ], [_βx₀, _βy₀, _βz₀]) #, _induced_current_coeffs)
    diodes = DIODES(_sensor_scale_factors, _azi_angles, _elev_angles)
    satellite_truth = SATELLITE(_J, magnetometer, diodes) 

    # Initialize modes/structs
    operation_mode = mag_cal    # First mode: Calibrate the magnetometer

    # Initialize histories to max lenth 
    x_hist = zeros(length(x0), _max_sim_length + 1)
    sensors_hist      = SENSORS( zeros(3, _max_sim_length), zeros(3, _max_sim_length), zeros(_num_diodes, _max_sim_length))
    ground_truth_hist = GROUND_TRUTH( zeros(_max_sim_length) , zeros(3, _max_sim_length), zeros(3, _max_sim_length))


    # Initialize satellite estimates (to average of distributions?)
    m = MAGNETOMETER(zeros(3), zeros(3), zeros(3))#, zeros(3,_num_diodes))
    d = DIODES(zeros(_num_diodes), zeros(_num_diodes), zeros(_num_diodes))
    _̂J = _J
    satellite_estimate = SATELLITE(_̂J, m, d) # Jhat = J for now

    # Record
    estimates_hist = Array{SATELLITE,1}(undef, _max_sim_length)
    for i = 1:_max_sim_length 
        new_sat = SATELLITE(_̂J, 
                    MAGNETOMETER(zeros(3), zeros(3), zeros(3)),# zeros(3,_num_diodes)), 
                    DIODES(zeros(_num_diodes), zeros(_num_diodes), zeros(_num_diodes)))
        estimates_hist[i] = new_sat
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
    ground_truth_hist.t_hist[1] = 0.0 # _epc
    


    @showprogress "Working" for i = 1:_max_sim_length
        state = x_hist[:, i]
        t = ground_truth_hist.t_hist[i] + (i - 1) * _dt

        # Add in SIM vs HARDWARE options
        truth, sensors = generate_measurements(sim, satellite_truth, albedo, state, t + _epc, CONSTANTS)

        
        #### In TRUTH, should they be UNIT or not (its mixed!)?

        # Could also use a Title_2_index[string(operation_mode)]  |  (e.g., "mag_cal" => index 2)
        estimators = sensors_to_estimators(sensors, truth)

        if true #i % 60 == 0
            satellite_estimate, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode])
        else
            estimator_finished = false 
        end
        ## AT SOME POINT I NEED TO USE MY MAG CALIB TO ADJUST THE MAG FIELD! DO IN MEKF x 2? or do it in the loop after calibration

        controllers = sensors_to_controllers(sensors)

        control_input, controller_finished = generate_command(controllers[operation_mode], state)

        new_state = rk4(satellite_truth, state, control_input, t + _epc, _dt)


        # Update mode - Currently ugly b/c cant iterate an @enum value... yet
        if (operation_mode == mag_cal) && (estimator_finished)
            operation_mode = finished
        elseif (operation_mode == detumble) && (controller_finished)
            operation_mode = diode_cal 
        elseif (operation_mode == diode_cal) && (estimator_finished)
            operation_mode = full_control 
        elseif (operation_mode == full_control) && (controller_finished)
            operation_mode = finished
        end

        # Update histories
        x_hist[:, i + 1] = new_state

        sensors_hist.magnetometer[:, i] = sensors.magnetometer 
        sensors_hist.sun[:, i] = sensors.sun 
        sensors_hist.diodes[:, i] = sensors.diodes 

        ground_truth_hist.t_hist[i] = t
        ground_truth_hist.Bᴵ_hist[:,i] = truth.Bᴵ_hist 
        ground_truth_hist.sᴵ_hist[:,i] = truth.sᴵ_hist 

        estimates_hist[i] = deepcopy(satellite_estimate)

        if operation_mode == finished
            x_hist = x_hist[:, 1:i]

            break
        end

    end

    temp = 1 # in case I need to return something extra for debugging
    return temp, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist
    
end

# UPDATE
function compute_diode_albedo(alb_mat, surf_normal, pos)
    """ Generates effect of Earth's albedo on one diode """
    return 1.0
end

# FIX
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH)
    # Assumes that IGRF IS PERFECT
    """ [MAG_CALIB; TRIVIAL; DIODE_CALIB; MEKF] """
    # MAG_CALIB: B_meas, B_pred, current 
    B_meas = sens.magnetometer
    B_pred = truth.Bᴵ_hist
    currents = sens.diodes 
    mag_calib = MAG_CALIB(B_meas, B_pred) #, currents)

    # TRIVIAL: Garbage 
    triv = TRIVIAL(0.0)

    # DIODE_CALIB: lots -> UPDATE!
    diode_calib = DIODE_CALIB(1.0, 1.0, 6)

    # MEKF: Lots, but less lots  -> Update
    mekf = MEKF(1.0)
    # return [mag_calib; triv; diode_calib; mekf]
    return [mag_calib; triv; triv; triv]
end

# UPDATE
function sensors_to_controllers(sens::SENSORS)
    """[Triv; detumbler; triv; mag]"""
    t = TRIVIAL(1.0)
    return [t; t; t; t]
end





temp, sat_true, sat_est, x_hist, ground_truth_hist, sensors_hist, estimates_hist = main();

display(plot(x_hist[1:3,:]', title = "Position"))

tru1 = plot(ground_truth_hist.Bᴵ_hist', title = "Bᴵ")
tru2 = plot(ground_truth_hist.sᴵ_hist', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1)))

sen1 = plot(sensors_hist.magnetometer', title = "Bᴮ")
sen2 = plot(sensors_hist.sun', title = "sᴮ")
sen3 = plot(sensors_hist.diodes', title = "Diodes")
display(plot(sen1, sen2, sen3, layout = (3,1)))


N = size(x_hist, 2) - 1
a, b, c = sat_true.magnetometer.scale_factors 
ρ, λ, ϕ = sat_true.magnetometer.non_ortho_angles
a_est, b_est, c_est = zeros(N), zeros(N), zeros(N)
ρ_est, λ_est, ϕ_est = zeros(N), zeros(N), zeros(N)
for i = 1:N
    a_est[i], b_est[i], c_est[i] = estimates_hist[i].magnetometer.scale_factors
    ρ_est[i], λ_est[i], ϕ_est[i] = estimates_hist[i].magnetometer.non_ortho_angles
end

ap = plot(a_est, title = "ã", label = "Estimate"); ap = hline!([a], linestyle = :dash, label = "Truth", ylims = [0, 2.0])
bp = plot(b_est, title = "b̃", label = "Estimate"); bp = hline!([b], linestyle = :dash, label = "Truth", ylims = [0, 2.0])
cp = plot(c_est, title = "c̃", label = "Estimate"); cp = hline!([c], linestyle = :dash, label = "Truth", ylims = [0, 2.0])
display(plot(ap, bp, cp, layout = (3, 1)))


ρp = plot(ρ_est, title = "ρ", label = "Estimate"); ap = hline!([ρ], linestyle = :dash, label = "Truth")
λp = plot(λ_est, title = "λ", label = "Estimate"); bp = hline!([λ], linestyle = :dash, label = "Truth")
ϕp = plot(ϕ_est, title = "ϕ", label = "Estimate"); cp = hline!([ϕ], linestyle = :dash, label = "Truth")
display(plot(ap, bp, cp, layout = (3, 1)))
