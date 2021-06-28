# Main file for mission simulator

# USING and INCLUDES that are used in ALL 
using Plots
using LinearAlgebra   # GENERIC/common ones
using SatelliteDynamics   # Allows for coordinate transforms, eclipse, sun position stuff, etc...
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


# Should I create a run() function that does initialization and then calls a loop/main() function
#   that takes in arguments? Then I can just adjust the arguments...?
function main()

    # Generate our satellite 
    magnetometer = MAGNETOMETER([_a, _b, _c], [_ρ, _λ, _ϕ], [_βx₀, _βy₀, _βz₀]) 
    diodes = DIODES(_sensor_scale_factors, _azi_angles, _elev_angles)
    satellite_truth = SATELLITE(_J, magnetometer, diodes) 

    # Initialize modes/structs
    operation_mode = mag_cal    # First mode: Calibrate the magnetometer

    # Initialize histories to max lenth 
    x_hist = zeros(length(x0), _max_sim_length + 1)
    # sensors_hist      = SENSORS( zeros(3, _max_sim_length), zeros(3, _max_sim_length), zeros(_num_diodes, _max_sim_length), zeros(3, _max_sim_length), zeros(3, _max_sim_length))
    # ground_truth_hist = GROUND_TRUTH( zeros(_max_sim_length) , zeros(3, _max_sim_length), zeros(3, _max_sim_length))


    # Initialize satellite estimates (to average of distributions?)
    m = MAGNETOMETER(ones(3), zeros(3), zeros(3)) # Assume no bias, no non-orthogonality, and a unit scale factor 
    m = initialize(m)

    d = DIODES(zeros(_num_diodes), zeros(_num_diodes), zeros(_num_diodes))
    d = initialize(satellite_truth, d) # Calls on 
    _̂J = _J
    satellite_estimate = SATELLITE(_̂J, m, d) # Jhat = J for now

    # Record
    estimates_hist = Array{SATELLITE,1}(undef, _max_sim_length)
    sensors_hist = Array{SENSORS, 1}(undef, _max_sim_length)
    ground_truth_hist = Array{GROUND_TRUTH, 1}(undef, _max_sim_length)
    for i = 1:_max_sim_length 
        new_sat = SATELLITE(_̂J, 
                    MAGNETOMETER(zeros(3), zeros(3), zeros(3)),# zeros(3,_num_diodes)), 
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
    # updated_data = MAG_CALIB(0.0, 0.0)


    operation_mode = diode_cal
    updated_data = DIODE_CALIB(albedo,
                                0.0, # NOT the same as x0
                                0.0, # Just empty stuff to be filled in 
                                0.0, # UPDATE with rᴵ
                                x0[11:13], 
                                0.0, # UPDATE with rᴮ
                                0.0, # UPDATE with current_meas 
                                0.0, # Update with W 
                                0.0, # Update with V 
                                SYSTEM._dt, 
                                SYSTEM._epc,
                                _num_diodes,
                                x0[1:3], # Position
                                true) # First pass

    satellite_estimate.magnetometer = satellite_truth.magnetometer

    temp = zeros(3, _max_sim_length) # Used for un-corrected mag field
    temp_cor = zeros(3, _max_sim_length)



    @showprogress "Working" for i = 1:_max_sim_length
        state = x_hist[:, i]
        t = ground_truth_hist[i].t_hist + (i - 1) * _dt

        # Add in SIM vs HARDWARE options
        #### In TRUTH, should they be UNIT or not (its mixed!)?
        truth, sensors = generate_measurements(sim, satellite_truth, albedo, state, t, CONSTANTS)

        temp[:, i] = sensors.magnetometer

        sensors.magnetometer = correct_mag_field(satellite_estimate, sensors.magnetometer) # Correct using our calibration estimate

        temp_cor[:, i] = sensors.magnetometer

        # Could also use a Title_2_index[string(operation_mode)]  |  (e.g., "mag_cal" => index 2)
        estimators = sensors_to_estimators(sensors, truth, updated_data)

        if true # i % 60 == 0
            satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode])
        else
            estimator_finished = false 
        end
        ## AT SOME POINT I NEED TO USE MY MAG CALIB TO ADJUST THE MAG FIELD! DO IN MEKF x 2? or do it in the loop after calibration

        controllers = sensors_to_controllers(sensors)

        control_input, controller_finished = generate_command(controllers[operation_mode], state)

        new_state = rk4(satellite_truth, state, control_input, t, _dt)


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

        # sensors_hist.magnetometer[:, i] = sensors.magnetometer 
        # sensors_hist.sun[:, i] = sensors.sun 
        # sensors_hist.diodes[:, i] = sensors.diodes 
        # sensors_hist.gyro[:, i] = sensors.gyro
        # sensors_hist.gps[:, i] = sensors.gps

        # ground_truth_hist.t_hist[i] = t
        # ground_truth_hist.Bᴵ_hist[:,i] = truth.Bᴵ_hist 
        # ground_truth_hist.sᴵ_hist[:,i] = truth.sᴵ_hist 

        sensors_hist[i] = SENSORS(sensors.magnetometer, sensors.sun, sensors.diodes, sensors.gyro, sensors.gps)
        ground_truth_hist[i] = GROUND_TRUTH(t, truth.Bᴵ_hist, truth.sᴵ_hist)
        estimates_hist[i] = deepcopy(satellite_estimate)

        if operation_mode == finished
            x_hist = x_hist[:, 1:i]

            break
        end

    end

    # temp = 1 # in case I need to return something extra for debugging
    return temp, temp_cor, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist
    
end


# FIX
function sensors_to_estimators(sens::SENSORS, truth::GROUND_TRUTH, data::MAG_CALIB)
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
    # diode_calib = DIODE_CALIB(1.0, 1.0, 6)

    # MEKF: Lots, but less lots  -> Update
    mekf = MEKF(1.0)
    # return [mag_calib; triv; diode_calib; mekf]
    return [mag_calib; triv; triv; triv]
end

# Do I need to return a set if I am overloading?
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
    # return [mag_calib; triv; diode_calib; mekf]
    return [mag_calib; triv; diode_calib; triv]
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

# UPDATE
function sensors_to_controllers(sens::SENSORS)
    """[Triv; detumbler; triv; mag]"""
    t = TRIVIAL(1.0)
    return [t; t; t; t]
end



temp, temp_cor, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist = main();

display(plot(x_hist[1:3, :]', title = "Position"))


# sen1 = plot(sensors_hist.magnetometer', title = "Bᴮ")
# sen2 = plot(sensors_hist.sun', title = "sᴮ")
# sen3 = plot(sensors_hist.diodes', title = "Diodes")
# sen4 = plot(sensors_hist.gyro', title = "Gyro")
# display(plot(sen1, sen2, sen3, sen4, layout = (4,1)))

N = size(x_hist, 2) - 1
c_est = zeros(6, N)
ϵ_est = zeros(6, N)
α_est = zeros(6, N)
Bᴵ = zeros(3, N)
sᴵ = zeros(3, N)
for i = 1:N
    c_est[:, i] = estimates_hist[i].diodes.calib_values 
    α_est[:, i] = estimates_hist[i].diodes.azi_angles
    ϵ_est[:, i] = estimates_hist[i].diodes.elev_angles
    Bᴵ[:, i] = ground_truth_hist[i].Bᴵ_hist
    sᴵ[:, i] = ground_truth_hist[i].sᴵ_hist
end

tru1 = plot(Bᴵ', title = "Bᴵ")
tru2 = plot(sᴵ', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1)))

c₁ = plot(c_est[1, :]); c₁ = hline!([satellite_truth.diodes.calib_values[1]], linestyle = :dash)
c₂ = plot(c_est[2, :]); c₂ = hline!([satellite_truth.diodes.calib_values[2]], linestyle = :dash)
c₃ = plot(c_est[3, :]); c₃ = hline!([satellite_truth.diodes.calib_values[3]], linestyle = :dash)
c₄ = plot(c_est[4, :]); c₄ = hline!([satellite_truth.diodes.calib_values[4]], linestyle = :dash)
c₅ = plot(c_est[5, :]); c₅ = hline!([satellite_truth.diodes.calib_values[5]], linestyle = :dash)
c₆ = plot(c_est[6, :]); c₆ = hline!([satellite_truth.diodes.calib_values[6]], linestyle = :dash)
display(plot(c₁, c₂, c₃, c₄, c₅, c₆, layout = (3, 2)))

a₁ = plot(α_est[1, :]); a₁ = hline!([satellite_truth.diodes.azi_angles[1]], linestyle = :dash)
a₂ = plot(α_est[2, :]); a₂ = hline!([satellite_truth.diodes.azi_angles[2]], linestyle = :dash)
a₃ = plot(α_est[3, :]); a₃ = hline!([satellite_truth.diodes.azi_angles[3]], linestyle = :dash)
a₄ = plot(α_est[4, :]); a₄ = hline!([satellite_truth.diodes.azi_angles[4]], linestyle = :dash)
a₅ = plot(α_est[5, :]); a₅ = hline!([satellite_truth.diodes.azi_angles[5]], linestyle = :dash)
a₆ = plot(α_est[6, :]); a₆ = hline!([satellite_truth.diodes.azi_angles[6]], linestyle = :dash)
display(plot(a₁, a₂, a₃, a₄, a₅, a₆, layout = (3, 2)))



# ##### MAGNETOMETER CALIBRATION
# temp, sat_true, sat_est, x_hist, ground_truth_hist, sensors_hist, estimates_hist = main();
# display(plot(x_hist[1:3,:]', title = "Position"))

# tru1 = plot(ground_truth_hist.Bᴵ_hist', title = "Bᴵ")
# tru2 = plot(ground_truth_hist.sᴵ_hist', title = "sᴵ")
# display(plot(tru1, tru2, layout = (2,1)))

# sen1 = plot(sensors_hist.magnetometer', title = "Bᴮ")
# sen2 = plot(sensors_hist.sun', title = "sᴮ")
# sen3 = plot(sensors_hist.diodes', title = "Diodes")
# sen4 = plot(sensors_hist.gyro', title = "Gyro")
# display(plot(sen1, sen2, sen3, sen4, layout = (4,1)))


# N = size(x_hist, 2) - 1
# a, b, c = sat_true.magnetometer.scale_factors 
# ρ, λ, ϕ = sat_true.magnetometer.non_ortho_angles
# a_est, b_est, c_est = zeros(N), zeros(N), zeros(N)
# ρ_est, λ_est, ϕ_est = zeros(N), zeros(N), zeros(N)
# for i = 1:N
#     a_est[i], b_est[i], c_est[i] = estimates_hist[i].magnetometer.scale_factors
#     ρ_est[i], λ_est[i], ϕ_est[i] = estimates_hist[i].magnetometer.non_ortho_angles
# end

# ap = plot(a_est, title = "ã", label = "Estimate"); ap = hline!([a], linestyle = :dash, label = "Truth", ylims = [0, 2.0])
# bp = plot(b_est, title = "b̃", label = "Estimate"); bp = hline!([b], linestyle = :dash, label = "Truth", ylims = [0, 2.0])
# cp = plot(c_est, title = "c̃", label = "Estimate"); cp = hline!([c], linestyle = :dash, label = "Truth", ylims = [0, 2.0])
# display(plot(ap, bp, cp, layout = (3, 1)))


# ρp = plot(ρ_est, title = "ρ", label = "Estimate"); ap = hline!([ρ], linestyle = :dash, label = "Truth")
# λp = plot(λ_est, title = "λ", label = "Estimate"); bp = hline!([λ], linestyle = :dash, label = "Truth")
# ϕp = plot(ϕ_est, title = "ϕ", label = "Estimate"); cp = hline!([ϕ], linestyle = :dash, label = "Truth")
# display(plot(ap, bp, cp, layout = (3, 1)))

# aᶠ = a_est[end]
# bᶠ = b_est[end]
# cᶠ = c_est[end]
# ρᶠ = ρ_est[end]
# λᶠ = λ_est[end]
# ϕᶠ = ϕ_est[end]
# println("a: |$a - $aᶠ| = ", a - aᶠ)
# println("b: |$b - $bᶠ| = ", b - bᶠ)
# println("c: |$c - $cᶠ| = ", c - cᶠ)
# println("ρ: |$ρ - $ρᶠ| = ", ρ - ρᶠ)
# println("λ: |$λ - $λᶠ| = ", λ - λᶠ)
# println("ϕ: |$ϕ - $ϕᶠ| = ", ϕ - ϕᶠ)

# T = [a          0               0;
#      b*sin(ρ)   b*cos(ρ)        0; 
#      c*sin(λ)   c*sin(ϕ)*cos(λ) c*cos(ϕ)*cos(λ)  ]
     

# Tᶠ = [aᶠ           0                  0;
#      bᶠ*sin(ρᶠ)   bᶠ*cos(ρᶠ)         0; 
#      cᶠ*sin(λᶠ)   cᶠ*sin(ϕᶠ)*cos(λᶠ) cᶠ*cos(ϕᶠ)*cos(λᶠ)  ]

# println(T)
# println(Tᶠ)