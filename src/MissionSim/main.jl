# [src/MissionSim/main.jl]

""" TO Do 
 - Make it not static 
 - Easy way to update sat by state (STATE and SAT_STATE)

 - Make this a module, return a dict from main, etc... (include state_machine)
"""

using Infiltrator, Test

using StaticArrays, SatelliteDynamics, EarthAlbedo 
using Distributions, LinearAlgebra, Plots, JLD2, Random
using ProgressMeter


include("mag_field.jl");     # Contains IGRF13 stuff 
include("quaternions.jl")

include("CustomStructs.jl");         using .CustomStructs
include("Simulator/Simulator.jl");   using .Simulator
include("Estimator/Estimator.jl");   using .Estimator 
include("Controller/Controller.jl"); using .Controller

@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, mekf, chill, finished)   
Base.to_index(om::Operation_mode) = Int(s)
include("state_machine.jl")



""" x₀ = initial state (eci is r, v); 
ℓ is max sim length """
# Add in commands for different starting modes, noise, etc... # NumOrbits, whether it is detumbled and bias_less
@info "Using FD in prediction"
@info "No noise, no albedo (?)"
function main(; t₀::Epoch = Epoch(2021, 1, 1), N = 6, dt = 0.2) 
    ##### INITIALIZE ##### 

    # Get initial state (enviroment)
    x₀, T_orbit = get_initial_state(; bias_less = false, detumbled = false)
    # x₀ = STATE(x₀.r, x₀.v, x₀.q, SVector{3, Float64}(0.0, deg2rad(0.5), 0.0), x₀.β); @info "Low spin!"

    #   Make sure that the initial environment state matches the true satellite state
    sat_state_truth = SAT_STATE(; q = x₀.q, β = x₀.β)
    sat_truth = SATELLITE(; sta = sat_state_truth)

    sat_est   = SATELLITE(; J = sat_truth.J, sta = SAT_STATE(; ideal = true), mag = MAGNETOMETER(; ideal = true), dio = DIODES(; ideal = true), cov = sat_truth.covariance)
    x, t = x₀, t₀
    ℓ = Int(round(4.0 * T_orbit / dt))  # MAX length of sim, in time steps

    alb = get_albedo(2);  # Set up REFL data

    # Allocate the history vectors
    truth, sensor, ecl, noise = generate_measurements(sat_truth, alb, x, t₀, dt; use_albedo = false)
    truths   = [truth   for _ = 1:ℓ]
    sensors  = [sensor  for _ = 1:ℓ]
    ecls     = zeros(ℓ)
    noises   = [noise   for _ = 1:ℓ]
    states   = [x       for _ = 1:ℓ]
    sat_ests = [sat_est for _ = 1:ℓ]
    modes    = [mag_cal for _ = 1:ℓ]

    ### TODO - Make this better. Starting conditions affect which part of satellite are ideal/match, ℓ, whether it is detumbled, etc...
    # Specify start conditions 
    flags = FLAGS(; init_detumble = false, mag_cal = false, dio_cal = false, final_detumble = false, in_sun = false)
    progress_bar = Progress(ℓ);

    ## Start with Detumble
    op_mode = detumble; data = nothing;  # Start with detumble

    ## Start with Mag Cal 
    # op_mode = detumble; flags = FLAGS(; init_detumble = true, mag_cal = false, dio_cal = false, final_detumble = false); data = nothing

    ## Start with Diode Cal 
    # op_mode = chill; flags.init_detumble = true; flags.magnetometer_calibrated = true; data = nothing;
    #     sat_est = SATELLITE(sat_est.J, sat_truth.magnetometer, sat_est.diodes, sat_est.state, sat_est.covariance)

    ## Start with Detumble (Round II)
    # op_mode = detumble; data = nothing; flags = FLAGS(; init_detumble = true, mag_cal = true, dio_cal = true, final_detumble = false)
    #     sat_est = SATELLITE(sat_est.J, sat_truth.magnetometer, sat_truth.diodes, sat_truth.state, sat_est.covariance)

    ## Start with MEKF 
    # op_mode = detumble; flags = FLAGS(; init_detumble = true, mag_cal = true, dio_cal = true, final_detumble = true); data = nothing
    #     sat_est = SATELLITE(sat_est.J, sat_truth.magnetometer, sat_truth.diodes, sat_est.state, sat_est.covariance);



    ### Call Loop 
    for i = 1:ℓ
        prev_mode = op_mode

        # Step
        sat_truth, sat_est, x, t, op_mode, data, truth, sensor, ecl, noise  = step(sat_truth, sat_est, alb, x, t, 
                                                                    dt, op_mode, flags, i, data, 
                                                                    progress_bar, T_orbit; use_albedo = false, σβ = 0.0)

        # Evaluate detumbling 
        (prev_mode == detumble) && (op_mode != detumble) && detumbler_report(states[1:i], sensors[1:i])

        # Evaluate performance of magnetometer calibration 
        (prev_mode ==   mag_cal) && (op_mode != mag_cal) && magnetometer_calibration_report(sat_truth, sat_est, sat_ests[1])

        # Evaluate performance of diode calibration 
        (prev_mode == diode_cal) && (op_mode != diode_cal) && (flags.diodes_calibrated) && diode_calibration_report(sat_truth, sat_ests[1:i-1]) 

        # Update histories
        truths[i]    = truth 
        sensors[i]   = sensor
        ecls[i]      = deepcopy(ecl) 
        noises[i]    = noise   # SAME PROBELM 
        states[i]    = x
        sat_ests[i]  = sat_est
        modes[i]     = op_mode # SAME PROBLEM 

        if op_mode == finished  # Trim the data and break 
            truths    = truths[1:i - 1]
            sensors   = sensors[1:i - 1]
            ecls      = ecls[1:i - 1] 
            noises    = noises[1:i - 1]
            states    = states[1:i - 1]
            sat_ests  = sat_ests[1:i - 1] 
            modes     = modes[1:i - 1] 

            @info "BREAKING EARLY!"
            break
        end
    end

    return sat_truth, sat_est, truths, sensors, ecls, noises, states, sat_ests, modes
end

function get_initial_state(; _Re = 6378136.3, detumbled = false, bias_less = false) 
    ecc = 0.0001717 + 0.00001 * randn()
    inc = 51.6426 + randn()
    Ω   = 178.1369 + randn()
    ω   = 174.7410 + randn()
    M   = 330.7918 + 100 + randn()   # +94/95 is just before sun, -40 is just before eclipse
    sma = (_Re + 421e3) / (1 + ecc)  # Apogee = semi_major * (1 + ecc)

    oe0 = [sma, ecc, inc, Ω, ω, M]   # Initial state, oscullating elements
    eci0 = sOSCtoCART(oe0, use_degrees = true) # Convert to Cartesean

    r₀ = SVector{3, Float64}(eci0[1:3])
    v₀ = SVector{3, Float64}(eci0[4:6])
    q₀ = randn(4);  q₀ = SVector{4, Float64}(q₀ / norm(q₀))
    ω₀ = (detumbled) ? SVector{3, Float64}(0.05 * randn(3)) : SVector{3, Float64}(0.5 * randn(3))
    β₀ = (bias_less) ? SVector{3, Float64}(0.0, 0.0, 0.0)  : SVector{3, Float64}(rand(Normal(0.0, deg2rad(2.0)), 3))
    
    T_orbit = orbit_period(oe0[1])
    x = STATE(r₀, v₀, q₀, ω₀, β₀)

    return x, T_orbit 
end

function get_albedo(scale = 1) 

    function load_refl(path = "data/refl.jld2", scale = 1)
        temp = load(path)
    
        refl = REFL( temp["data"][1:scale:end, 1:scale:end], temp["type"], temp["start_time"], temp["stop_time"])
    
        return refl
    end
    lat_step = 1.0 * scale
    lon_step = 1.25 * scale

    refl = load_refl("data/refl.jld2", scale)  
    cell_centers_ecef = get_albedo_cell_centers(lat_step, lon_step) 
    return Simulator.ALBEDO(refl, cell_centers_ecef)
end;

function get_albedo_cell_centers(lat_step = 1, lon_step = 1.25)
    """
        Returns the cell centers for the grid covering the surface of the Earth in Cartesian ECEF, to be used in later estimations of Earth's albedo,
            by looping through each cell's LLA coordinate and converting to ECEF 

        Arguments:
        - lat_step: (Optional) The step size (in degrees) to take across the latitude domain. Defaults to 1*        | Scalar 
        - lon_step: (Optional) The step size (in degrees) to take across the longitude domain. Defaults to 1.25*    | Scalar

        Returns:
        - cells_ecef: Matrix containing [x,y,z] coordinate for each latitude, longitude point.
                        Of form [lat, lon, [x,y,z]]                                                                 | [num_lat x num_lon x 3]
    """
    alt = 0.0 # Assume all cells are on surface of earth
    num_lat = Int(round((180 - lat_step) / lat_step) + 1)
    num_lon = Int(round((360 - lon_step) / lon_step) + 1)

    lon_offset = lon_step + (360 - lon_step) / 2   # Centers at 0 (longitude: [1.25, 360] => [-179.375, 179.375])
    lat_offset = lat_step + (180 - lat_step) / 2   # Centers at 0 (latitude:  [1.00, 180] => [-89.5, 89.5])

    cells_ecef = zeros(num_lat, num_lon, 3) # Lat, Lon, [x,y,z]
    for lat = 1:num_lat 
        for lon = 1:num_lon
            geod = [(lon * lon_step - lon_offset), (lat * lat_step - lat_offset), alt]
            ecef = sGEODtoECEF(geod, use_degrees = true)

            cells_ecef[Int(lat), Int(lon), :] = ecef
        end
    end

    return cells_ecef 
end;





# Leave in main 
function magnetometer_calibration_report(sat_true, sat_est, sat_init) 

    aᶠ, bᶠ, cᶠ = round.(sat_est.magnetometer.scale_factors, sigdigits = 3)
    ρᶠ, λᶠ, ϕᶠ = round.(sat_est.magnetometer.non_ortho_angles, sigdigits = 3)
    βxᶠ, βyᶠ, βzᶠ = round.(sat_est.magnetometer.bias, sigdigits = 3)

    a, b, c = round.(sat_true.magnetometer.scale_factors, sigdigits = 3)
    ρ, λ, ϕ = round.(sat_true.magnetometer.non_ortho_angles, sigdigits = 3)
    βx, βy, βz = round.(sat_true.magnetometer.bias, sigdigits = 3)

    a₀, b₀, c₀ = round.(sat_init.magnetometer.scale_factors, sigdigits = 3)
    ρ₀, λ₀, ϕ₀ = round.(sat_init.magnetometer.non_ortho_angles, sigdigits = 3)
    βx₀, βy₀, βz₀ = round.(sat_init.magnetometer.bias, sigdigits = 3)

    println("__________________________________________________________________________")
    println("___PARAM___|___Truth____|__Final Guess__|__Init Guess__|__Improved?__")
    println("     a     |   $a\t|    $aᶠ  \t|    $a₀       | ", abs(a - aᶠ) < abs(a - a₀) ? "    Yes!" : "   No!")
    println("     b     |   $b\t|    $bᶠ  \t|    $b₀       | ", abs(b - bᶠ) < abs(b - b₀) ? "    Yes!" : "   No!")
    println("     c     |   $c\t|    $cᶠ  \t|    $c₀       | ", abs(c - cᶠ) < abs(c - c₀) ? "    Yes!" : "   No!")
    println("     ρ°    |   $ρ\t|    $ρᶠ  \t|    $ρ₀       | ", abs(ρ - ρᶠ) < abs(ρ - ρ₀) ? "    Yes!" : "   No!")
    println("     λ°    |   $λ\t|    $λᶠ  \t|    $λ₀       | ", abs(λ - λᶠ) < abs(λ - λ₀) ? "    Yes!" : "   No!")
    println("     ϕ°    |   $ϕ\t|    $ϕᶠ  \t|    $ϕ₀       | ", abs(ϕ - ϕᶠ) < abs(ϕ - ϕ₀) ? "    Yes!" : "   No!")
    println("     βx    |   $βx\t|    $βxᶠ \t|    $βx₀       | ", abs(βx - βxᶠ) < abs(βx - βx₀) ? "    Yes!" : "   No!")
    println("     βy    |   $βy\t|    $βyᶠ \t|    $βy₀       | ", abs(βy - βyᶠ) < abs(βy - βy₀) ? "    Yes!" : "   No!")
    println("     βz    |   $βz\t|    $βzᶠ \t|    $βz₀       | ", abs(βz - βzᶠ) < abs(βz - βz₀) ? "    Yes!" : "   No!")
    println("__________________________________________________________________________")

    return nothing
end;

# Error in sun vec estimation before and after?
function diode_calibration_report(sat_true::SATELLITE, est_hist::Vector{SATELLITE{6, T}}) where {T} 
    N = size(est_hist, 1)
    C_est, α_est, ϵ_est = zeros(6, N), zeros(6, N), zeros(6, N)

    for i = 1:N
        C_est[:, i] = est_hist[i].diodes.calib_values
        α_est[:, i] = rad2deg.(est_hist[i].diodes.azi_angles)
        ϵ_est[:, i] = rad2deg.(est_hist[i].diodes.elev_angles)
    end
    
    C₀, α₀, ϵ₀ = est_hist[1].diodes.calib_values, rad2deg.(est_hist[1].diodes.azi_angles), rad2deg.(est_hist[1].diodes.elev_angles)
    Cf, αf, ϵf = est_hist[N].diodes.calib_values, rad2deg.(est_hist[N].diodes.azi_angles), rad2deg.(est_hist[N].diodes.elev_angles)
    C, α, ϵ = sat_true.diodes.calib_values, rad2deg.(sat_true.diodes.azi_angles), rad2deg.(sat_true.diodes.elev_angles)
    
    println("_____________________________________DIODES______________________________________")
    println("____DIODE___|______Truth (C,α,ϵ)_______|______Final Guess______|______Init Guess_______|_____Improved?___")
    for i = 1:6 
        print("     $i     |    $(round(C[i], digits = 2)), $(round(α[i], digits = 2)), $(round(ϵ[i], digits = 2))  \t|")
        print("   $(round(Cf[i], digits = 2)), $(round(αf[i], digits = 2)), $(round(ϵf[i], digits = 2))  \t|")
        print("   $(round(C₀[i], digits = 2)), $(round(α₀[i], digits = 2)), $(round(ϵ₀[i], digits = 2))  \t|")
        print( (abs(C₀[i] - C[i]) ≤ abs(Cf[i] - C[i])) ? " No!, " : " Yes!, ")
        print( (abs(α₀[i] - α[i]) ≤ abs(αf[i] - α[i])) ? " No!, " : " Yes!, ")
        println( (abs(ϵ₀[i] - ϵ[i]) ≤ abs(ϵf[i] - ϵ[i])) ? " No! " : " Yes! ")
    end
    println("__________________________________________________________________________________")
    

    ##### PLOT #####
    C_off, ϵ_off, α_off = 0.2, 5.0, 5.0
    Cps, ϵps, αps = [], [], []
    for i = 1:6
        plot(C_est[i, :], title = "Scale Factor (C)", label = false); 
            hline!([C₀[i]], ls = :dot, label = false); Cp = hline!([C[i]], ls = :dash, label = false, ylim = [C[i] - C_off, C[i] + C_off])
        
        plot(ϵ_est[i, :], title = "Elevation Angle (ϵ)", label = false); 
            hline!([ϵ₀[i]], ls = :dot, label = false); ϵp = hline!([ϵ[i]], ls = :dash, label = false, ylim = [ϵ[i] - ϵ_off, ϵ[i] + ϵ_off])
        
        plot(α_est[i, :], title = "Azimuth Angle (α)", label = false); 
            hline!([α₀[i]], ls = :dot, label = false); αp = hline!([α[i]], ls = :dash, label = false, ylim = [α[i] - α_off, α[i] + α_off])
    

        push!(Cps, Cp)
        push!(ϵps, ϵp)
        push!(αps, αp)
    end

    # Adjust the labels 
    Cps[2].series_list[1][:label] = "Estimate"; Cps[2].series_list[2][:label] = "Initial Guess"; Cps[2].series_list[3][:label] = "Truth"
    ϵps[2].series_list[1][:label] = "Estimate"; ϵps[2].series_list[2][:label] = "Initial Guess"; ϵps[2].series_list[3][:label] = "Truth"
    αps[2].series_list[1][:label] = "Estimate"; αps[2].series_list[2][:label] = "Initial Guess"; αps[2].series_list[3][:label] = "Truth"

    display(plot(Cps..., layout = (3, 2)))
    display(plot(ϵps..., layout = (3, 2)))
    display(plot(αps..., layout = (3, 2)))
end;

function mekf_report(states::Vector{STATE{T}}, est_hist::Vector{SATELLITE{6, T}}) where {T}
    N = size(states, 1)
    qs = [states[i].q for i = 1:N]; qs = reduce(hcat, qs)'; 
    βs = [states[i].β for i = 1:N]; βs = reduce(hcat, βs)';
    q̂s = [est_hist[i].state.q for i = 1:N]; q̂s = reduce(hcat, q̂s)';
    β̂s = [est_hist[i].state.β for i = 1:N]; β̂s = reduce(hcat, β̂s)';

    qErrs = [norm(cayley_map(qs[i, :], q̂s[i, :])) for i = 1:N];
    βErrs = [norm(βs[i, :] - β̂s[i, :]) for i = 1:N];

    plot( qs, title = "MEKF Report: q", c = [:red :orange :blue :green]);
    plot!(q̂s, c = [:red :orange :blue :green], ls = :dash, label = false);
    display( plot!(qErrs, label = false, c = :black, ls = :dot, ylim = (-1.5, 1.5)) )

    plot( βs, title = "MEKF Report: β", c = [:red :blue :green]);
    plot!(β̂s, c = [:red :blue :green], ls = :dash, label = false);
    display( plot!(βErrs, label = false, c = :black, ls = :dot) )
 
    return qs, q̂s, βs, β̂s
end;

function detumbler_report(states, sensors; τ₁ = deg2rad(25), τ₂ = deg2rad(10) )
    N = size(states, 1)
    ωs = [states[i].ω for i = 1:N]; ωs = reduce(hcat, ωs)'
    ω̂s = [sensors[i].gyro for i = 1:N]; ω̂s = reduce(hcat, ω̂s)'
    nω = [norm(ωs[i, :]) for i = 1:N];
    nω̂ = [norm(ω̂s[i, :]) for i = 1:N];

    plot(ωs, c = [:red :blue :green], label = ["ωx" "ωy" "ωz"]); hline!([τ₁], c = :gray, ls = :dot, label = "τ₀"); 
        hline!([τ₂], c = :gray, ls = :dot, label = "τf"); 
        display(plot!(nω, ls = :dash, c = :black, label = "Mag", xlabel = "Step", ylabel = "Ang Velocity (rad/s)", title = "True Angular Velocity"))

    plot(ω̂s, c = [:red :blue :green], label = ["ωx" "ωy" "ωz"]); hline!([τ₁], c = :gray, ls = :dot, label = "τ₀"); 
        hline!([τ₂], c = :gray, ls = :dot, label = "τf"); 
        display(plot!(nω̂, ls = :dash, c = :black, label = "Mag", xlabel = "Step", ylabel = "Ang Velocity (rad/s)", title = "Est Angular Velocity"))

    return ωs, ω̂s
end;



# Random.seed!(1001)
sat_truth, sat_est, truths, sensors, ecls, noises, states, sat_ests, op_modes = main(); # Adjust arguments to include initial state
# mekf_report(states, sat_ests)
# diode_calibration_report(sat_truth, sat_ests)

println("Done");
