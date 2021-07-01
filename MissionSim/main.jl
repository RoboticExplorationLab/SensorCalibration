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


# Should I create a run() function that does initialization and then calls a loop/main() function  -> YES
#   that takes in arguments? Then I can just adjust the arguments...?  -> MAYBE
function main()

    println("Initializing Values")

    # Generate our satellite 
    magnetometer = MAGNETOMETER([_a, _b, _c], [_ρ, _λ, _ϕ], [_βx₀, _βy₀, _βz₀]) 
    diodes = DIODES(_sensor_scale_factors, _azi_angles, _elev_angles)
    satellite_truth = SATELLITE(_J, magnetometer, diodes) 

    # Initialize modes/structs
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
    updated_data = initialize(albedo, x0, SYSTEM)

    satellite_estimate.magnetometer = satellite_truth.magnetometer  ###############

    operation_mode = detumble

    temp = zeros(3, _max_sim_length) # Used for un-corrected mag field
    temp_cor = zeros(3, _max_sim_length)
    count = 0



    @showprogress "Working" for i = 1:_max_sim_length
        state = x_hist[:, i]
        t = ground_truth_hist[i].t_hist + (i - 1) * _dt

        # Add in SIM vs HARDWARE options
        #### In TRUTH, should they be UNIT or not (its mixed!)?
        truth, sensors = generate_measurements(sim, satellite_truth, albedo, state, t, CONSTANTS, _dt)

        temp[:, i] = sensors.magnetometer

        sensors.magnetometer = correct_mag_field(satellite_estimate, sensors.magnetometer) # Correct using our calibration estimate

        temp_cor[:, i] = sensors.magnetometer

        # Could also use a Title_2_index[string(operation_mode)]  |  (e.g., "mag_cal" => index 2)
        estimators = sensors_to_estimators(sensors, truth, updated_data)

        old_estimate = deepcopy(satellite_estimate)
        if true #i % 60 == 0
            satellite_estimate, updated_data, estimator_finished = estimate_vals(satellite_estimate, estimators[operation_mode])
        else
            estimator_finished = false 
        end

        if operation_mode == diode_cal 
            estimator_finished = check_if_finished(old_estimate, satellite_estimate, updated_data)
            if count < 4000
                estimator_finished = false 
            end
        end
        # estimator_finished = false 

        controllers = sensors_to_controllers(sensors, _dt)

        control_input, controller_finished = generate_command(controllers[operation_mode], state)
        if operation_mode == detumble
            controller_finished = check_if_finished(x_hist, i)
        end

        new_state = rk4(satellite_truth, state, control_input, t, _dt)

        # Update mode - Currently ugly b/c cant iterate an @enum value... yet
        if (operation_mode == mag_cal) && (estimator_finished)
            operation_mode = detumble
            updated_data = TRIVIAL(0.0)
        elseif (operation_mode == detumble) && (controller_finished)
            println("\nFinished Detumbling")
            operation_mode = diode_cal 
            count = 0
            updated_data = initialize(albedo, new_state, SYSTEM)
        elseif (operation_mode == diode_cal) && (estimator_finished)
            println("\nFinished Calibrating Diodes")
            operation_mode = finished 
        elseif (operation_mode == full_control) && (controller_finished)
            operation_mode = finished
        end

        # Update histories
        x_hist[:, i + 1] = new_state
        x_hist[7:10, i + 1] = x_hist[7:10, i + 1] / norm(x_hist[7:10, i + 1])

        sensors_hist[i] = SENSORS(sensors.magnetometer, sensors.sun, sensors.diodes, sensors.gyro, sensors.gps)
        ground_truth_hist[i] = GROUND_TRUTH(t, truth.Bᴵ_hist, truth.sᴵ_hist)
        estimates_hist[i] = deepcopy(satellite_estimate)

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
    return [mag_calib; triv; triv; triv]
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
    
    return [triv; triv; diode_calib; triv]
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





temp, temp_cor, satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist = main();

display(plot(x_hist[1:3, :]', title = "Position"))


DETUMBLER
N = size(x_hist, 2) - 1
ŵ = zeros(3, N)
mag_w = zeros(N)

for i = 1:N 
    ŵ[:, i] = rad2deg.(sensors_hist[i].gyro)
    mag_w[i] = rad2deg(norm(x_hist[11:13, i]))
end

w = plot(rad2deg.(x_hist[11:13, :])', color = [:red :blue :green], label = ["wₓ" "wy" "wz"], title = "Ang Vel (rad/sec)")
w = plot!(ŵ', linestyle = :dash, color = [:red :blue :green], label = ["ŵₓ" "ŵy" "ŵz"])
w = plot!(mag_w, linestyle = :dot, color = :black, label = "|w_true|")
display(w)

bias = x_hist[14:16, :]

display(plot(bias', title = "Bias (random walk)"))



##### STUFF FOR PLOTTING DIODE CALIBRATION
N = size(x_hist, 2) - 1
c_est = zeros(6, N)
ϵ_est = zeros(6, N)
α_est = zeros(6, N)
Bᴵ = zeros(3, N)
sᴵ = zeros(3, N)
for i = 1:N
    c_est[:, i] = estimates_hist[i].diodes.calib_values
    α_est[:, i] = rad2deg.(estimates_hist[i].diodes.azi_angles)
    ϵ_est[:, i] = rad2deg.(estimates_hist[i].diodes.elev_angles)
    Bᴵ[:, i] = ground_truth_hist[i].Bᴵ_hist
    sᴵ[:, i] = ground_truth_hist[i].sᴵ_hist
end

tru1 = plot(Bᴵ', title = "Bᴵ")
tru2 = plot(sᴵ', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1)))

c₁ = plot(c_est[1, :], ylabel = "Scale Factor Estimate", xlabel = "Time", title = "Scale Factors"); c₁ = hline!([satellite_truth.diodes.calib_values[1]], ylim = [satellite_truth.diodes.calib_values[1] - 0.25, satellite_truth.diodes.calib_values[1] + 0.25], linestyle = :dash); c₁ = hline!([c_est[1,1]], linestyle = :dot)
c₂ = plot(c_est[2, :]); c₂ = hline!([satellite_truth.diodes.calib_values[2]], ylim = [satellite_truth.diodes.calib_values[2] - 0.25, satellite_truth.diodes.calib_values[2] + 0.25], linestyle = :dash); c₂ = hline!([c_est[2,1]], linestyle = :dot)
c₃ = plot(c_est[3, :]); c₃ = hline!([satellite_truth.diodes.calib_values[3]], ylim = [satellite_truth.diodes.calib_values[3] - 0.25, satellite_truth.diodes.calib_values[3] + 0.25],linestyle = :dash); c₃ = hline!([c_est[3,1]], linestyle = :dot)
c₄ = plot(c_est[4, :]); c₄ = hline!([satellite_truth.diodes.calib_values[4]], ylim = [satellite_truth.diodes.calib_values[4] - 0.25, satellite_truth.diodes.calib_values[4] + 0.25],linestyle = :dash); c₄ = hline!([c_est[4,1]], linestyle = :dot)
c₅ = plot(c_est[5, :]); c₅ = hline!([satellite_truth.diodes.calib_values[5]], ylim = [satellite_truth.diodes.calib_values[5] - 0.25, satellite_truth.diodes.calib_values[5] + 0.25],linestyle = :dash); c₅ = hline!([c_est[5,1]], linestyle = :dot)
c₆ = plot(c_est[6, :]); c₆ = hline!([satellite_truth.diodes.calib_values[6]], ylim = [satellite_truth.diodes.calib_values[6] - 0.25, satellite_truth.diodes.calib_values[6] + 0.25],linestyle = :dash); c₆ = hline!([c_est[6,1]], linestyle = :dot)
display(plot(c₁, c₂, c₃, c₄, c₅, c₆, layout = (3, 2)))

offset = 7.5
a₁ = plot(α_est[1, :], ylabel = "Azi Estimate", xlabel = "Time", title = "Azimuth Angles"); a₁ = hline!([rad2deg(satellite_truth.diodes.azi_angles[1])], ylim = [rad2deg(satellite_truth.diodes.azi_angles[1]) - offset, rad2deg(satellite_truth.diodes.azi_angles[1]) + offset], linestyle = :dash); a₁ = hline!([α_est[1,1]], linestyle = :dot)
a₂ = plot(α_est[2, :]); a₂ = hline!([rad2deg(satellite_truth.diodes.azi_angles[2])], ylim = [rad2deg(satellite_truth.diodes.azi_angles[2]) - offset, rad2deg(satellite_truth.diodes.azi_angles[2]) + offset], linestyle = :dash); a₂ = hline!([α_est[2,1]], linestyle = :dot)
a₃ = plot(α_est[3, :]); a₃ = hline!([rad2deg(satellite_truth.diodes.azi_angles[3])], ylim = [rad2deg(satellite_truth.diodes.azi_angles[3]) - offset, rad2deg(satellite_truth.diodes.azi_angles[3]) + offset],  linestyle = :dash); a₃ = hline!([α_est[3,1]], linestyle = :dot)
a₄ = plot(α_est[4, :]); a₄ = hline!([rad2deg(satellite_truth.diodes.azi_angles[4])], ylim = [rad2deg(satellite_truth.diodes.azi_angles[4]) - offset, rad2deg(satellite_truth.diodes.azi_angles[4]) + offset],  linestyle = :dash); a₄ = hline!([α_est[4,1]], linestyle = :dot)
a₅ = plot(α_est[5, :]); a₅ = hline!([rad2deg(satellite_truth.diodes.azi_angles[5])], ylim = [rad2deg(satellite_truth.diodes.azi_angles[5]) - offset, rad2deg(satellite_truth.diodes.azi_angles[5]) + offset],  linestyle = :dash); a₅ = hline!([α_est[5,1]], linestyle = :dot)
a₆ = plot(α_est[6, :]); a₆ = hline!([rad2deg(satellite_truth.diodes.azi_angles[6])], ylim = [rad2deg(satellite_truth.diodes.azi_angles[6]) - offset, rad2deg(satellite_truth.diodes.azi_angles[6]) + offset],  linestyle = :dash); a₆ = hline!([α_est[6,1]], linestyle = :dot)
display(plot(a₁, a₂, a₃, a₄, a₅, a₆, layout = (3, 2)))

e₁ = plot(ϵ_est[1, :], ylabel = "Elev Estimate", xlabel = "Time", title = "Elevation Angles"); a₁ = hline!([rad2deg(satellite_truth.diodes.elev_angles[1])], ylim = [rad2deg(satellite_truth.diodes.elev_angles[1]) - offset, rad2deg(satellite_truth.diodes.elev_angles[1]) + offset],  linestyle = :dash); e₁ = hline!([ϵ_est[1,1]], linestyle = :dot)
e₂ = plot(ϵ_est[2, :]); a₂ = hline!([rad2deg(satellite_truth.diodes.elev_angles[2])], ylim = [rad2deg(satellite_truth.diodes.elev_angles[2]) - offset, rad2deg(satellite_truth.diodes.elev_angles[2]) + offset], linestyle = :dash); e₂ = hline!([ϵ_est[2,1]], linestyle = :dot)
e₃ = plot(ϵ_est[3, :]); a₃ = hline!([rad2deg(satellite_truth.diodes.elev_angles[3])], ylim = [rad2deg(satellite_truth.diodes.elev_angles[3]) - offset, rad2deg(satellite_truth.diodes.elev_angles[3]) + offset], linestyle = :dash); e₃ = hline!([ϵ_est[3,1]], linestyle = :dot)
e₄ = plot(ϵ_est[4, :]); a₄ = hline!([rad2deg(satellite_truth.diodes.elev_angles[4])], ylim = [rad2deg(satellite_truth.diodes.elev_angles[4]) - offset, rad2deg(satellite_truth.diodes.elev_angles[4]) + offset], linestyle = :dash); e₄ = hline!([ϵ_est[4,1]], linestyle = :dot)
e₅ = plot(ϵ_est[5, :]); a₅ = hline!([rad2deg(satellite_truth.diodes.elev_angles[5])], ylim = [rad2deg(satellite_truth.diodes.elev_angles[5]) - offset, rad2deg(satellite_truth.diodes.elev_angles[5]) + offset], linestyle = :dash); e₅ = hline!([ϵ_est[5,1]], linestyle = :dot)
e₆ = plot(ϵ_est[6, :]); a₆ = hline!([rad2deg(satellite_truth.diodes.elev_angles[6])], ylim = [rad2deg(satellite_truth.diodes.elev_angles[6]) - offset, rad2deg(satellite_truth.diodes.elev_angles[6]) + offset], linestyle = :dash); e₆ = hline!([ϵ_est[6,1]], linestyle = :dot)
display(plot(e₁, e₂, e₃, e₄, e₅, e₆, layout = (3, 2)))



# ##### MAGNETOMETER CALIBRATION
# N = size(x_hist, 2) - 1
# Bᴵ = zeros(3, N)
# sᴵ = zeros(3, N)

# diodes = zeros(6, N)
# gyro = zeros(3, N)
# Bᴮ = zeros(3, N)
# sᴮ = zeros(3, N)
# for i = 1:N
#     Bᴮ[:, i] = sensors_hist[i].magnetometer
#     sᴮ[:, i] = sensors_hist[i].sun
#     gyro[:, i] = sensors_hist[i].gyro 
#     diodes[:, i] = sensors_hist[i].diodes
#     Bᴵ[:, i] = ground_truth_hist[i].Bᴵ_hist
#     sᴵ[:, i] = ground_truth_hist[i].sᴵ_hist
# end
# tru1 = plot(Bᴵ', title = "Bᴵ")
# tru2 = plot(sᴵ', title = "sᴵ")
# display(plot(tru1, tru2, layout = (2,1)))

# sen1 = plot(Bᴮ', title = "Bᴮ")
# sen2 = plot(sᴮ', title = "sᴮ")
# sen3 = plot(diodes', title = "Diodes")
# sen4 = plot(gyro', title = "Gyro")
# display(plot(sen1, sen2, sen3, sen4, layout = (4,1)))


# a, b, c = satellite_truth.magnetometer.scale_factors 
# ρ, λ, ϕ = satellite_truth.magnetometer.non_ortho_angles
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