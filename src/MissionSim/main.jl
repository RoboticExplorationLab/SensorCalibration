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

using Logging; logger = SimpleLogger(stdout, Logging.Debug);


include("mag_field.jl");     # Contains IGRF13 stuff 
include("quaternions.jl")

include("CustomStructs.jl");         using .CustomStructs
include("Simulator/Simulator.jl");   using .Simulator
include("Estimator/Estimator.jl");   using .Estimator 
include("Controller/Controller.jl"); using .Controller

@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, mekf, chill, finished)   
Base.to_index(om::Operation_mode) = Int(s)
RecipesBase.plot(m::Vector{Operation_mode}) =  plot(Int.(m), yticks = ([1:6;], ["Mag", "Det", "Dio", "MEKF", "Chi", "Fin"]))

include("state_machine.jl")
include("reports.jl")



""" x₀ = initial state (eci is r, v); 
ℓ is max sim length """
# Add in verbose for showing plots or not (used for monte carlo)
function main(; t₀::Epoch = Epoch(2021, 1, 1), N = 6, dt = 0.2, verbose = true, initial_state::Operation_mode = detumble, 
                num_orbits = 3, albedo_ds = 2, use_albedo = true, kwargs...) 


    ##### INITIALIZE ##### 

    ### Set up desired starting state
    if initial_state == detumble 
        """ Default initial state. Sat est initializes to all ideals"""
        x₀, T_orbit = get_initial_state(; detumbled = false)
        flags = FLAGS(; init_detumble = false, mag_cal = false, dio_cal = false, final_detumble = false, in_sun = false)
        sat_truth = SATELLITE(; sta =  SAT_STATE(; q = x₀.q, β = x₀.β))
        sat_est   = SATELLITE(; J = sat_truth.J, sta = SAT_STATE(; ideal = true), mag = MAGNETOMETER(; ideal = true), dio = DIODES(; ideal = true))
        op_mode   = detumble

    elseif initial_state == mag_cal 
        """ Starts with a low angular velocity """
        x₀, T_orbit = get_initial_state(; detumbled = true)
        flags = FLAGS(; init_detumble = true)
        sat_truth = SATELLITE(; sta =  SAT_STATE(; q = x₀.q, β = x₀.β))
        sat_est   = SATELLITE(; J = sat_truth.J, sta = SAT_STATE(; ideal = true), mag = MAGNETOMETER(; ideal = true), dio = DIODES(; ideal = true))
        op_mode   = detumble  # Will switch over after first iteration

    elseif initial_state == diode_cal
        """ Starts with a low angular velocity and correct magnetometer parameters """
        x₀, T_orbit = get_initial_state(; detumbled = true)
        flags = FLAGS(; init_detumble = true, mag_cal = true)
        sat_truth = SATELLITE(; sta =  SAT_STATE(; q = x₀.q, β = x₀.β))
        sat_est   = SATELLITE(; J = sat_truth.J, mag = sat_truth.magnetometer, dio = DIODES(; ideal = true), sta = SAT_STATE(; ideal = true))
        op_mode = chill;      # Will switch over after first iteration

    elseif initial_state == mekf 
        """ Starts with a low angular velocity and correct magnetometer/diode parameters """
        x₀, T_orbit = get_initial_state(; detumbled = true)
        flags = FLAGS(; init_detumble = true, mag_cal = true, dio_cal = true, final_detumble = true)
        sat_truth = SATELLITE(; sta =  SAT_STATE(; q = x₀.q, β = x₀.β))
        sat_est   = SATELLITE(; J = sat_truth.J, mag = sat_truth.magnetometer, dio = sat_truth.diodes, sta = SAT_STATE(; ideal = true))
        op_mode = chill;      # Will switch over after first iteration

    else
        @error "Invalid initial state! Cannot use $initial_state"
    end

    ℓ = Int(round(num_orbits * T_orbit / dt))   # MAX length of sim, in time steps
    alb = get_albedo(albedo_ds);                # Set up REFL data
    progress_bar = Progress(ℓ);                 # Set up progress bar


    ### Allocate the history vectors
    truth, sensor, ecl, noise = generate_measurements(sat_truth, alb, x₀, t₀, dt; use_albedo = false) # Not used, just for sizing
    truths   = [truth   for _ = 1:ℓ]
    sensors  = [sensor  for _ = 1:ℓ]
    ecls     = zeros(ℓ)
    noises   = [noise   for _ = 1:ℓ]
    states   = [x₀      for _ = 1:ℓ]
    sat_ests = [sat_est for _ = 1:ℓ]
    modes    = [mag_cal for _ = 1:ℓ]


    ℓ = Int(round(num_orbits * T_orbit / dt))  # MAX length of sim, in time steps
    alb = get_albedo(albedo_ds);  # Set up REFL data
    progress_bar = Progress(ℓ);

    ### Call Loop 
    data = nothing    # Will be updated in state machine
    x, t = x₀, t₀
    for i = 1:ℓ
        prev_mode = op_mode
            
        ########################## add in noise args here
        # Step
        sat_truth, sat_est, x, t, op_mode, data, truth, sensor, ecl, noise  = step(sat_truth, sat_est, alb, x, t, 
                                                                        dt, op_mode, flags, i, progress_bar, T_orbit, data; 
                                                                        use_albedo = use_albedo, kwargs...)

        if verbose
            # Evaluate detumbling 
            (prev_mode == detumble) && (op_mode != detumble) && detumbler_report(states[1:i - 1], sensors[1:i - 1])

            # Evaluate performance of magnetometer calibration 
            (prev_mode ==   mag_cal) && (op_mode != mag_cal) && magnetometer_calibration_report(sat_truth, sat_est, sat_ests[1])

            # Evaluate performance of diode calibration 
            (prev_mode == diode_cal) && (op_mode != diode_cal) && (flags.diodes_calibrated) && diode_calibration_report(sat_truth, sat_ests[1:i-1]) 

            # Evaluate performance of MEKF
            (prev_mode == mekf) && (op_mode != mekf) && (flags.diodes_calibrated) && mekf_report(states[1:i-1], sat_ests[1:i-1]) 
        end

        # Update histories
        truths[i]    = truth 
        sensors[i]   = sensor
        ecls[i]      = deepcopy(ecl) 
        noises[i]    = noise   
        states[i]    = x
        sat_ests[i]  = sat_est
        modes[i]     = op_mode 

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

    results = (sat_truth = sat_truth, 
            sat_ests  = sat_ests,
            truths    = truths,
            sensors   = sensors, 
            ecls      = ecls,
            states    = states, 
            modes     = modes, 
            noises    = noises)
    
    return results
end

function get_initial_state(; _Re = 6378136.3, detumbled = false, bias_less = false) 
    ecc = 0.0001717 + 0.001 * randn()
    inc = 51.6426 + 5 * randn()
    Ω   = 178.1369 + 5 * randn()
    ω   = 174.7410 + 5 * randn()
    M   = 330.7918 + 50 * randn()   # +94/95 is just before sun, -40 is just before eclipse
    sma = (_Re + 421e3 + 1000 * randn()) / (1 + ecc)  # Apogee = semi_major * (1 + ecc)

    oe0 = [sma, ecc, inc, Ω, ω, M]   # Initial state, oscullating elements
    eci0 = sOSCtoCART(oe0, use_degrees = true) # Convert to Cartesean

    r₀ = SVector{3, Float64}(eci0[1:3])
    v₀ = SVector{3, Float64}(eci0[4:6])
    q₀ = randn(4);  q₀ = SVector{4, Float64}(q₀ / norm(q₀))

    # If it is too low, it is bad for calibration, so we don't want zero mean
    ω₀ = (detumbled) ? rand(Normal(0.06, 0.02)) : rand(Normal(0.35, 0.15), 3)
    ω₀ = SVector{3, Float64}(ω₀ .* sign.(randn(3)))
    
    # ω₀ = (detumbled) ? SVector{3, Float64}(0.07 * randn(3)) : SVector{3, Float64}(0.4 * randn(3))
    β₀ = (bias_less) ? SVector{3, Float64}(0.0, 0.0, 0.0)  : SVector{3, Float64}(rand(Normal(0.0, deg2rad(2)), 3)) # Initial guess can be a bit off
    
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






# nd = [norm(sensors[i].diodes) for i = 1:size(sensors, 1)]
# ndc = [norm(sensors[i].diodes ./ sat_ests[i].diodes.calib_values) for i = 1:size(sensors, 1)]

# for i = 1:6 
#     d = [sensors[j].diodes[i] for j = 1000:3000]
#     plot!(d)
# end

# with_logger(logger) do 
#     global temp = main(; num_orbits = 0.5, initial_state = diode_cal)
# end



# @info "No Noise!"; results = main(; σβ = 0.0, σB = 0.0, σ_gyro = 0.0, σr = 0.0, σ_current = 0.0); 
# sat_truth, sat_est, truths, sensors, ecls, noises, states, sat_ests, op_modes 

# display(plot(results[:states]))
# display(plot(results[:sensors]))
# diode_calibration_report(results)
# mekf_report(results)
