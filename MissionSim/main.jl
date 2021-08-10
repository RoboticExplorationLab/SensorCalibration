# Main file for mission simulator

using Plots
using LinearAlgebra, SatelliteDynamics        
using Random, Distributions
using ProgressMeter
using EarthAlbedo
using MAT, JLD2
using StaticArrays

using Infiltrator
# Random.seed!(12344) # 1234  | Alt: 4321

include("mag_field.jl")            # Contains IGRF13 stuff
include("rotationFunctions.jl")    # Contains general functions for working with quaternions

include("CustomStructs.jl");         using .CustomStructs   # 
include("system_config_file.jl"); 

include("Estimator/Estimator.jl");   using .Estimator
include("Simulator/Simulator.jl");   using .Simulator
include("Controller/Controller.jl"); using .Controller

include("state_machine.jl")

__init__() # Generates PyCall commands

@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, full_control, mekf, finished, chill)
Base.to_index(s::Operation_mode) = Int(s)

function report_on_magnetometer(truth::SATELLITE, est::SATELLITE, init::SATELLITE, sens_hist, truth_hist)
    aᶠ, bᶠ, cᶠ = round.(est.magnetometer.scale_factors, sigdigits = 3)
    ρᶠ, λᶠ, ϕᶠ = round.(est.magnetometer.non_ortho_angles, sigdigits = 3)
    βxᶠ, βyᶠ, βzᶠ = round.(est.magnetometer.bias, sigdigits = 3)

    a, b, c = round.(truth.magnetometer.scale_factors, sigdigits = 3)
    ρ, λ, ϕ = round.(truth.magnetometer.non_ortho_angles, sigdigits = 3)
    βx, βy, βz = round.(truth.magnetometer.bias, sigdigits = 3)

    a₀, b₀, c₀ = round.(init.magnetometer.scale_factors, sigdigits = 3)
    ρ₀, λ₀, ϕ₀ = round.(init.magnetometer.non_ortho_angles, sigdigits = 3)
    βx₀, βy₀, βz₀ = round.(init.magnetometer.bias, sigdigits = 3)

    println("__________________________________________________________________________")
    println("___PARAM___|___Truth____|__Final Guess__|__Init Guess__|__Improved?__")
    println("     a     |   $a\t|    $aᶠ  \t|    $a₀       | ", abs(a - aᶠ) < abs(a - a₀) ? "better!" : "worse!")
    println("     b     |   $b\t|    $bᶠ  \t|    $b₀       | ", abs(b - bᶠ) < abs(b - b₀) ? "better!" : "worse!")
    println("     c     |   $c\t|    $cᶠ  \t|    $c₀       | ", abs(c - cᶠ) < abs(c - c₀) ? "better!" : "worse!")
    println("     ρ°    |   $ρ\t|    $ρᶠ  \t|    $ρ₀       | ", abs(ρ - ρᶠ) < abs(ρ - ρ₀) ? "better!" : "worse!")
    println("     λ°    |   $λ\t|    $λᶠ  \t|    $λ₀       | ", abs(λ - λᶠ) < abs(λ - λ₀) ? "better!" : "worse!")
    println("     ϕ°    |   $ϕ\t|    $ϕᶠ  \t|    $ϕ₀       | ", abs(ϕ - ϕᶠ) < abs(ϕ - ϕ₀) ? "better!" : "worse!")
    println("     βx    |   $βx\t|    $βxᶠ \t|    $βx₀       | ", abs(βx - βxᶠ) < abs(βx - βx₀) ? "better!" : "worse!")
    println("     βy    |   $βy\t|    $βyᶠ \t|    $βy₀       | ", abs(βy - βyᶠ) < abs(βy - βy₀) ? "better!" : "worse!")
    println("     βz    |   $βz\t|    $βzᶠ \t|    $βz₀       | ", abs(βz - βzᶠ) < abs(βz - βz₀) ? "better!" : "worse!")
    println("__________________________________________________________________________")

    N = size(sens_hist, 1)
    B̃ᴮ = zeros(3, N)
    Bᴮ = zeros(3, N)
    Bᴮ_cor = zeros(3, N)

    a, b, c = est.magnetometer.scale_factors 
    ρ, λ, ϕ = est.magnetometer.non_ortho_angles 
    β = est.magnetometer.bias

    T̂ = [a          0               0;
         b*sin(ρ)   b*cos(ρ)        0; 
         c*sin(λ)   c*sin(ϕ)*cos(λ) c*cos(ϕ)*cos(λ)  ]

    for i = 1:N 
        B̃ᴮ[:, i] = sens_hist[i].magnetometer 
        Bᴮ[:, i] = truth_hist[i].Bᴮ_hist
        Bᴮ_cor[:, i] = T̂^(-1) * (B̃ᴮ[:, i] - β)
    end


    println("Plotting!")
    og = plot(B̃ᴮ', title = "Original", label = ["x" "y" "z"], color = [:red :blue :green])
    og = plot!(Bᴮ', label =  false, color = [:red :blue :green], linestyle = :dash)

    err_og = B̃ᴮ - Bᴮ 
    err_og_mag  = [norm(B̃ᴮ[:, i] - Bᴮ[:, i]) for i = 1:N]
    erro = plot(err_og', title = "Errors", label = ["x" "y" "z"], color = [:red :blue :green])
        erro = plot!(err_og_mag, linewidth = 1, linestyle = :dash, color = :black)

    display(plot(og, erro, layout = (2,1)))
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
    w = plot!(mag_w, linestyle = :dot, color = :black, label = "|w_true|", legend = false)
    display(w)

    bias = x[14:16, :]

    display(plot(bias', title = "Bias (random walk)", legend = false))
    
end
function report_on_diodes(truth::SATELLITE, est_hist::Array{SATELLITE,1})
    N = size(est_hist, 1)
    c_est, α_est, ϵ_est = zeros(6, N), zeros(6, N), zeros(6, N)

    for i = 1:N
        c_est[:, i] = est_hist[i].diodes.calib_values
        α_est[:, i] = rad2deg.(est_hist[i].diodes.azi_angles)
        ϵ_est[:, i] = rad2deg.(est_hist[i].diodes.elev_angles)
    end

    ϵ₀ = rad2deg.([(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)])
    α₀ = rad2deg.([0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi])
    c₀ = ones(6)

    cᶠ, ϵᶠ, αᶠ = c_est[:, end], ϵ_est[:, end], α_est[:, end]
    c, ϵ, α = truth.diodes.calib_values, rad2deg.(truth.diodes.elev_angles), rad2deg.(truth.diodes.azi_angles)

    println("\n----- DIODE REPORT -----")
    for i = 1:_num_diodes
        print("Diode $i - ")
        print( abs(c₀[i] - c[i]) ≤ abs(cᶠ[i] - c[i]) ? "C: Worse! " : "C: Better!"); print("  (", round(c[i] ,sigdigits = 3), " & ", round(cᶠ[i], sigdigits = 3), " vs ", round(c₀[i], sigdigits = 3) ," (init))    \t| ")
        print( abs(ϵ₀[i] - ϵ[i]) ≤ abs(ϵᶠ[i] - ϵ[i]) ? "E: Worse! " : "E: Better!"); print("  (", round(ϵ[i] ,sigdigits = 3), " & ", round(ϵᶠ[i], sigdigits = 3), " vs ", round(ϵ₀[i], sigdigits = 3) ," (init)) \t| ")
        print( abs(α₀[i] - α[i]) ≤ abs(αᶠ[i] - α[i]) ? "A: Worse! " : "A: Better!"); print("  (", round(α[i] ,sigdigits = 3), " & ", round(αᶠ[i], sigdigits = 3), " vs ", round(α₀[i], sigdigits = 3) ," (init))\n")
    end
    println("------------------------")

    c_off = 0.2
    e_off = 2.0
    a_off = 2.0
    cp, ap, ep = Array{Plots.Plot{Plots.GRBackend}, 1}(undef, _num_diodes), Array{Plots.Plot{Plots.GRBackend}, 1}(undef, _num_diodes), Array{Plots.Plot{Plots.GRBackend}, 1}(undef, _num_diodes)
    for i = 1:_num_diodes
        cp[i] = plot(c_est[i, :], title = "Scale Factor (C)", label = false); cp[i] = hline!([c₀[i]], linestyle = :dot, label = false); cp[i] = hline!([c[i]], ylim = [c[i] - c_off, c[i] + c_off], linestyle = :dash, label = false)
        ep[i] = plot(ϵ_est[i, :], title = "Elevation Angle (ϵ)", label = false); ep[i] = hline!([ϵ₀[i]], linestyle = :dot, label = false); ep[i] = hline!([ϵ[i]], ylim = [ϵ[i] - e_off, ϵ[i] + e_off], linestyle = :dash, label = false); 
        ap[i] = plot(α_est[i, :], title = "Azimuth Angle (α)", label = false); ap[i] = hline!([α₀[i]], linestyle = :dot, label = false); ap[i] = hline!([α[i]], ylim = [α[i] - a_off, α[i] + a_off], linestyle = :dash, label = false); 
    end
    cp[2].series_list[1][:label] = "Estimate"; cp[2].series_list[2][:label] = "Initial Guess"; cp[2].series_list[3][:label] = "Truth"
    ep[2].series_list[1][:label] = "Estimate"; ep[2].series_list[2][:label] = "Initial Guess"; ep[2].series_list[3][:label] = "Truth"
    ap[2].series_list[1][:label] = "Estimate"; ap[2].series_list[2][:label] = "Initial Guess"; ap[2].series_list[3][:label] = "Truth"

    display(plot(cp..., layout = (3, 2))); savefig(string(run_folder, "scale_factors.png"))
    display(plot(ep..., layout = (3, 2))); savefig(string(run_folder, "elev.png"))
    display(plot(ap..., layout = (3, 2))); savefig(string(run_folder, "azi.png"))
end
function report_on_mekf(mekf_hist, state) 
    N = size(mekf_hist, 1)
    poses, biases = zeros(4, N), zeros(3, N)

    for i = 1:N
        poses[:, i] = state[7:10, i]
        biases[:,i] = state[14:16,i]
    end
end
function mat_from_vec(a)  # FROM KEVIN / Attitude.jl
    "Turn a vector of vectors into a matrix"


    rows = length(a[1])
    columns = length(a)
    A = zeros(rows,columns)

    for i = 1:columns
        A[:,i] = a[i]
    end

    return A
end


run_folder = "garbage/"  
@info "Running in $run_folder"

function generate_satellite()
    # GENERATE THE TRUE SATELLITE
    magnetometer = MAGNETOMETER([_a, _b, _c], [_ρ, _λ, _ϕ], [_βx₀, _βy₀, _βz₀]) 
    diodes = DIODES(_sensor_scale_factors, _azi_angles, _elev_angles)
    initial_state = [q0[:]; β0[:]]
    initial_covariance = NaN * ones((3 + 3 + 3 * _num_diodes))
    initial_covariance = diagm(initial_covariance)
    satellite_truth = SATELLITE(_J, magnetometer, diodes, initial_state, initial_covariance) 

    # GENERATE THE INITIAL ESTIMATE
    m = initialize_magnetometer()
    d = initialize_diodes(_num_diodes) 
    _̂J = _J
    state_est = [0 0 0 1 0 0 0]  # Assume no bias and unit rotation
    satellite_estimate = SATELLITE(_̂J, m, d, state_est, initial_covariance) # Jhat = J for now

    return satellite_truth, satellite_estimate
end

function generate_histories(sat_est)
    x_hist = zeros(length(x0), _max_sim_length + 1)  
    estimates_hist = Array{SATELLITE,1}(undef, _max_sim_length)
    sensors_hist = Array{SENSORS, 1}(undef, _max_sim_length)
    ground_truth_hist = Array{GROUND_TRUTH, 1}(undef, _max_sim_length)
    noise_hist = Array{NOISE, 1}(undef, _max_sim_length)

    # estimates_hist = [deepcopy(sat_est) for i = 1:_max_sim_length]
    # sensors_hist = [SENSORS(zeros(3), zeros(_num_diodes), zeros(3), zeros(3)) for i = 1:_max_sim_length]


    for i = 1:_max_sim_length 
        new_sat = deepcopy(sat_est) 
        estimates_hist[i] = new_sat

        new_sens = SENSORS(zeros(3), zeros(_num_diodes), zeros(3), zeros(3))
        sensors_hist[i] = new_sens 

        new_truth = GROUND_TRUTH(_epc, zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3))
        ground_truth_hist[i] = new_truth

        new_noise = NOISE(zeros(6), zeros(3), zeros(3), zeros(3,3), zeros(3,3))
        noise_hist[i] = new_noise
    end

    return x_hist, estimates_hist, sensors_hist, ground_truth_hist, noise_hist
end

function main()
    satellite_truth, satellite_estimate = generate_satellite()

    x_hist, estimates_hist, sensors_hist, ground_truth_hist, noise_hist = generate_histories(satellite_estimate)

    # Generate ALBEDO data 
    #   refl_dict = matread("../../Earth_Albedo_Model/Processed_Data/tomsdata2005/2005/ga050101-051231.mat")
    refl_dict = matread("refl.mat")     # Same as ^ but saved in local directory
    refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])

    cell_centers_ecef = get_albedo_cell_centers()
    albedo = ALBEDO(refl, cell_centers_ecef)

    sim = SIM(1.0)
    x_hist[:, 1] = x0

    __init__()  # Inits python functions for PyCall

    progress_bar = Progress(_max_sim_length) 

    # Temporary histories for debugging
    ecl_hist = zeros(_max_sim_length)
    mode_hist = zeros(_max_sim_length)

    operation_mode = mag_cal    
    # updated_data = MAG_CALIB(0.0, 0.0)
    updated_data = initialize(albedo, x0, SYSTEM)
    satellite_estimate.magnetometer = satellite_truth.magnetometer 
    # satellite_estimate.diodes = satellite_truth.dides
    flags = FLAGS(false, true, false, false, false) # in_sun, mag_cal, diodes_cal, detumbling, calibrating 

    count = 0
    try 
        for i = 1:_max_sim_length
            state = x_hist[:, i]
            t = ground_truth_hist[i].t_hist + (i - 1) * _dt

            truth, sensors, ecl, noise = generate_measurements(sim, satellite_truth, albedo, state, t, CONSTANTS, _dt)

            old_op = operation_mode
            
            operation_mode, controller, estimator, flags = update_operation_mode(flags, sensors, SYSTEM, albedo, updated_data, t, satellite_estimate)
            if old_op != operation_mode
                println("Switched from $old_op to $operation_mode at $i ")
                if (old_op == mag_cal) & (i > 3000)
                    report_on_magnetometer(satellite_truth, satellite_estimate, estimates_hist[1], sensors_hist[1:i], ground_truth_hist[1:i])
                end
            end

            if operation_mode == mag_cal 
                if i % _ds_rate == 0 
                    satellite_estimate, updated_data = estimate_vals(satellite_estimate, estimator, (i > Int(round(2 * orbit_period(oe0[1]))/_dt)))  # i > Int(round(2...)) rather than "false" ?
                end 
            else
                satellite_estimate, updated_data = estimate_vals(satellite_estimate, estimator)
            end

            control_input = generate_command(controller)

            new_state = rk4(satellite_truth, state, control_input, t, _dt)
            new_state[7:10] /= norm(new_state[7:10]) # Normalize quaternions

            x_hist[:, i + 1] = new_state

            ### Update histories
            if flags.magnetometer_calibrated
                sensors.magnetometer = correct_mag_field(satellite_estimate, sensors.magnetometer)
            end
            
            sensors_hist[i] = SENSORS(sensors.magnetometer, sensors.diodes, sensors.gyro, sensors.gps) 
            ground_truth_hist[i] = GROUND_TRUTH(t, truth.Bᴵ_hist, truth.sᴵ_hist, truth.ŝᴮ_hist, truth.Bᴮ_hist)
            estimates_hist[i] = deepcopy(satellite_estimate)
            ecl_hist[i] = ecl
            noise_hist[i] = noise
            mode_hist[i] = Int(operation_mode)

            if operation_mode == finished
                x_hist = x_hist[:, 1:(i+1)]   # Save through most recent 
                # save_global_variables(satellite_truth, estimates_hist[1])  # Used to test magnetometer downsampling
                break
            end

            # Display stuff
            change = 0
            ecl = (norm(truth.ŝᴮ_hist) < 0.02) ? true : false
            temp = norm(satellite_estimate.covariance)

            disp_elev = round.(rad2deg.(satellite_estimate.diodes.elev_angles), sigdigits = 3)
            disp_azi  = round.(rad2deg.(satellite_estimate.diodes.azi_angles ), sigdigits = 3)
            disp_cur  = round.(sensors.diodes, sigdigits = 3)
            if operation_mode == diode_cal
                notes = string("Ecl: $ecl \t| ||COV||: $temp   \t| Diode Currents: $disp_cur") #ϵ: $disp_elev   \t| α: $disp_azi")
            elseif operation_mode == detumble 
                notes = string("Ecl: $ecl \t| ||ω||: ",  rad2deg(norm(new_state[11:13])), "\t| ||w̃||: ", rad2deg(norm(sensors.gyro)))
            else
                notes = string("Ecl: $ecl \t|  ", norm(truth.ŝᴮ_hist))
            end
            ProgressMeter.next!(progress_bar; showvalues = [(:Mode, operation_mode), (:Iteration, i), (:Notes, notes)])

            count += 1
        end
    catch e 
        if e isa InterruptException 
            println("Function terminated by user")
            x_hist = x_hist[:, 1:(count)]
            return satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist, ecl_hist, noise_hist
        else
          @show e     
          rethrow(e)
          x_hist = x_hist[:, 1:(count)]
          return satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist, ecl_hist, noise_hist
        end
    end

    return satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist, ecl_hist, noise_hist, mode_hist
end

satellite_truth, satellite_estimate, x_hist, ground_truth_hist, sensors_hist, estimates_hist, ecl_hist, noise_hist, mode_hist = main();

println("Finished")

N = size(x_hist, 2) - 1

report_on_magnetometer(satellite_truth, satellite_estimate, estimates_hist[1], sensors_hist[1:N], ground_truth_hist[1:N])
report_on_detumbler(x_hist, sensors_hist[1:N])
report_on_diodes(satellite_truth, estimates_hist[1:N])

Bᴵ = mat_from_vec([ground_truth_hist[i].Bᴵ_hist for i = 1:N])
sᴵ = mat_from_vec([ground_truth_hist[i].sᴵ_hist for i = 1:N])
ŝᴮ = mat_from_vec([ground_truth_hist[i].ŝᴮ_hist for i = 1:N])
Bᴮ_truth = mat_from_vec([ground_truth_hist[i].Bᴮ_hist for i = 1:N])
t  = [ground_truth_hist[i].t_hist - _epc for i = 1:N]
ground_truth_hist = nothing; GC.gc()

Bᴮ =       mat_from_vec([sensors_hist[i].magnetometer for i = 1:N])
currents = mat_from_vec([sensors_hist[i].diodes for i = 1:N])
gyro =     mat_from_vec([sensors_hist[i].gyro   for i = 1:N])
pos =      mat_from_vec([sensors_hist[i].gps    for i = 1:N])
sᴮ_ests  = mat_from_vec([estimate_sun_vector(sensors_hist[i], estimates_hist[i]) for i = 1:N])
sensors_hist = nothing; GC.gc()

β_ests = mat_from_vec([estimates_hist[i].state[5:7] for i = 1:N])
q_ests = mat_from_vec([estimates_hist[i].state[1:4] for i = 1:N])
c_est = mat_from_vec([estimates_hist[i].diodes.calib_values for i = 1:N])
α_est = mat_from_vec([rad2deg.(estimates_hist[i].diodes.azi_angles)  for i = 1:N])
ϵ_est = mat_from_vec([rad2deg.(estimates_hist[i].diodes.elev_angles) for i = 1:N])
C_cov = [norm(estimates_hist[i].covariance[7:12,  7:12]) for i = 1:N]
α_cov = [norm(estimates_hist[i].covariance[13:18,13:18]) for i = 1:N]
ϵ_cov = [norm(estimates_hist[i].covariance[19:24,19:24]) for i = 1:N]
total_cov = [norm(estimates_hist[i].covariance[7:24, 7:24]) for i = 1:N] # Total covariance for calibration states
estimates_hist = nothing; GC.gc()

diode_noise = mat_from_vec([noise_hist[i].diodes for i = 1:N])
gyro_noise  = mat_from_vec([noise_hist[i].gyro for i = 1:N])
gps_noise   = mat_from_vec([noise_hist[i].gps for i = 1:N])
noise_hist = nothing; GC.gc()

sᴵ_mag = [norm(sᴵ[:, i]) for i = 1:N]; sᴵ_mag /= maximum(sᴵ_mag);
sᴮ_mag = [norm(ŝᴮ[:, i]) for i = 1:N]
 
C, ϵ, α = satellite_truth.diodes.calib_values, rad2deg.(satellite_truth.diodes.elev_angles), rad2deg.(satellite_truth.diodes.azi_angles)

a, b, c = satellite_truth.magnetometer.scale_factors 
ρ, λ, ϕ = satellite_truth.magnetometer.non_ortho_angles



# @save string(run_folder, "data.jld2")


# PLOT
tru1 = plot(Bᴵ', title = "Bᴵ")
tru2 = plot(sᴵ', title = "sᴵ")
display(plot(tru1, tru2, layout = (2,1))); savefig(string(run_folder, "inertial_truth.png"))

sen1 = plot(Bᴮ', title = "Bᴮ")
sen2 = plot(ŝᴮ', title = "sᴮ")
sen3 = plot(currents', title = "Diodes")
sen4 = plot(gyro', title = "Gyro")
display(plot(sen1, sen2, sen3, sen4, layout = (4,1), legend = false))


# qp1 = plot(mekf_hist[1:4, 1:N]', title = "Attitude Estimate"); qp2 = plot(x_hist[7:10,:]', title = "Truth", label = false); 
# qp3 = plot( (x_hist[7:10,2:(end)] - mekf_hist[1:4, 1:N])', ylim = (-1.1, 1.1), title = "Dif1", label = false); qp4 = plot( (x_hist[7:10,2:(end)] + mekf_hist[1:4, 1:N])', ylim = (-1.1, 1.1), title = "Dif2", label = false); 
# display(plot(qp1, qp2, qp3, qp4, layout = (4,1)))

# bp1 = plot(mekf_hist[5:7, 1:N]',   title = "Bias Estimate"  ); bp2 = plot(x_hist[14:16,:]', title = "Truth"); bp3 = plot( (x_hist[14:16,2:(end)] - mekf_hist[5:7, 1:N])', title = "Dif", ylim = [-0.04, 0.04]); 
# display(plot(bp1, bp2, bp3, layout = (3,1))); savefig(string(run_folder, "bias.png"))


d1 = plot(currents[1,:], title = "Diode Currents")
d2 = plot(currents[2,:], title = "Diode Currents")
d3 = plot(currents[3,:], title = "Diode Currents")
d4 = plot(currents[4,:], title = "Diode Currents")
d5 = plot(currents[5,:], title = "Diode Currents")
d6 = plot(currents[6,:], title = "Diode Currents")
display(plot(d1, d2, d3, d4, d5, d6, layout = (3,2))); savefig(string(run_folder, "currents.png"))
display(plot(diode_noise', title = "Diode Noises"))

cc = plot(C_cov, title = "C Cov")
ac = plot(α_cov, title = "α Cov")
ec = plot(ϵ_cov, title = "ϵ Cov")
totc = plot(total_cov, title = "Total")
display(plot(cc, ac, ec, totc, layout = (4, 1), title = "Covariances")); savefig(string(run_folder, "cov.png"))

a = plot(ŝᴮ', title = "Truth"); b= plot(sᴮ_ests', title = "sᴮ Est"); c = plot((ŝᴮ - sᴮ_ests)', title = "Dif"); 
display(plot(b, a, c, layout = (3,1))); savefig(string(run_folder, "body_sun_vec.png"))

a = plot(sᴵ_mag, title = "sᴵ Mag", label = "sᴵ"); b = plot(sᴮ_mag, title = "sᴮ Mag", label = "sᴮ"); c = plot(ecl_hist[1:N], title = "Eclipse", label = "Eclipse")    
display(plot(a, b, c, layout = (3, 1))); 

ŵ = gyro - β_ests
w = x_hist[11:13,  1:N]
a = plot(ŵ', title = "ω Est"); b = plot(w', title = "Truth", label = false); c = plot((w - ŵ)', title = "Diff", label = false); d = plot(gyro_noise', title = "Noise", label = false)
display(plot(a, b, d, c, layout = (4,1))); savefig(string(run_folder, "gyro.png"))

pos_hat = pos 
pos = x_hist[1:3, 1:N]
a = plot(pos_hat', title = "Pos Est"); b = plot(pos', title = "Truth", label = false); c = plot((pos - pos_hat)', title = "Diff", label = false); d = plot(gps_noise', title = "Noise", label = false)
display(plot(a, b, d, c, layout = (4,1))); savefig(string(run_folder, "gps.png"))

a = plot(Bᴮ', title = "Bᴮ Est"); b = plot(Bᴮ_truth', title = "Truth", label = false); c = plot((Bᴮ_truth - Bᴮ)', title = "Diff", label = false); #d = plot(mag_noise', title = "Noise", label = false)
display(plot(a,b, c, layout = (3,1))); savefig(string(run_folder, "body_mag.png"))



# # Sensors: Truth + Mag, Est + Mag, Noise + mag, percentage of truth : noise 
# # Gyro 
# ŵ = gyro - β_ests; ŵ_mag = [norm(ŵ[:, i]) for i = 1:N]
# w = x_hist[11:13,  1:N]; w_mag = [norm(w[:, i]) for i = 1:N]
# gyro_noise_mag = [norm(gyro_noise[:, i]) for i = 1:N]
# gyro_error = (w - ŵ); gyro_error_mag = [norm(gyro_error[:, i]) for i = 1:N]
# a = plot(ŵ', title = "ω Est"); a = plot!(w_mag, color = :black)
# b = plot(w', title = "Truth", label = false); b = plot!(ŵ_mag, label = false, color = :black)
# c = plot(gyro_noise', title = "Noise", label = false); c = plot!(gyro_noise_mag, label = false, color = :black)
# d = plot(gyro_error', title = "Error", label = false); d = plot!(gyro_error_mag, label = false, color = :black)
# e = plot(w_mag / gyro_noise_mag)
# display(plot(a, b, c, d, e, layout = (5, 1)))

# # GPS 
# pos_mag = [norm(pos[:, i]) for i = 1:N]
# pos_t = x_hist[11:13,  1:N]; pos_t_mag = [norm(pos_t[:, i]) for i = 1:N]
# gps_noise_mag = [norm(gps_noise[:, i]) for i = 1:N]
# pos_error = (pos_t - pos); pos_error_mag = [norm(pos_error[:, i]) for i = 1:N]
# a = plot(pos', title = "Pos Est"); a = plot!(pos_mag, color = :black)
# b = plot(pos_t', title = "Truth", label = false); b = plot!(pos_t_mag, label = false, color = :black)
# c = plot(gps_noise', title = "Noise", label = false); c = plot!(gps_noise_mag, label = false, color = :black)
# d = plot(pos_error', title = "Error", label = false); d = plot!(pos_error_mag, label = false, color = :black)
# e = plot(pos_t_mag / gps_noise_mag)
# display(plot(a, b, c, d, e, layout = (5, 1)))

# est = pos;                  est_mag = [norm(est[:, i]) for i = 1:N]
# tru = x_hist[11:13, 1:N];   tru_mag = [norm(tru[:, i]) for i = 1:N]
# noise = gps_noise;          noise_mag = [norm(noise[:, i]) for i = 1:N] 
# err = (tru - est);          error_mag = [norm(err[:, i]) for i = 1:N]
# a = plot(est', title = "Mag Est");                  a = plot!(est_mag, color = :black)
# b = plot(tru', title = "Truth", label = false);     b = plot!(tru_mag, color = :black, label = false)
# c = plot(noise', title = "Noise", label = false);   c = plot!(noise_mag, color = :black, label = false)
# d = plot(err', title = "Error", label = false);     d = plot!(error_mag, color=:black, label = false)
# e = plot(tru_mag / noise_mag)
# display(plot(a, b, c, d, e layout = (5, 1)))


# # Magnetometer 
# est = Bᴮ; est_mag = [norm(est[:, i]) for i = 1:N]
# tru = Bᴮ_truth; tru_mag = [norm(tru[:, i]) for i = 1:N]
# # noise = 
# err = (tru - est); error_mag = [norm(err[:, i]) for i = 1:N]
# a = plot(est', title = "Mag Est");                  a = plot!(est_mag, color = :black)
# b = plot(tru', title = "Truth", label = false);     b = plot!(tru_mag, color = :black, label = false)
# d = plot(err', title = "Error", label = false);     d = plot!(error_mag, color=:black, label = false)
# display(plot(a, b, d, layout = (3, 1)))
