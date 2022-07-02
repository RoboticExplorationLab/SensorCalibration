"""
    TEST FILE for comparing Simulator.jl to simulator.py
"""

using Infiltrator, PyCall, Test, MAT 
using SatelliteDynamics, LinearAlgebra, EarthAlbedo 
using Random, Distributions

include("../mag_field.jl")            # Contains IGRF13 stuff
include("../rotationFunctions.jl");   # Contains general functions for working with attitude and quaternions
include("../CustomStructs.jl");         using .CustomStructs 
include("../system_config_file.jl"); 
include("Simulator.jl");   using .Simulator
include("../Estimator/Estimator.jl"); using .Estimator

function __init_python__()
    # Allows us to run python scripts with PyCall
    py"""
    import sys 
    sys.path.insert(0, "./") # Allows us to include nearby python scripts

    import numpy as np 

    from earthAlbedo import EarthAlbedo 
    from satelliteDynamics import *
    from simulator import Simulator
    """
end

function generate_satellite()
    # Generates a pair of satellites to use during testing

    # GENERATE THE TRUE SATELLITE
    magnetometer = MAGNETOMETER([_a, _b, _c], [_ρ, _λ, _ϕ], [_βx₀, _βy₀, _βz₀]) 
    diodes = DIODES(_sensor_scale_factors, _azi_angles, _elev_angles)
    initial_state = [q0[:]; β0[:]]
    cov = NaN * ones((3 + 3 + 3 * _num_diodes))
    initial_covariance = diagm(cov)

    ######### "Reset_covariance"
    cov[1:6] = (10*pi/180) * ones(6)
    initial_covariance = diagm(cov)
    #########

    satellite_truth = SATELLITE(_J, magnetometer, diodes, initial_state, initial_covariance) 

    # GENERATE THE INITIAL ESTIMATE
    m = initialize_magnetometer()
    d = initialize_diodes(_num_diodes) 
    _̂J = _J
    state_est = [0 0 0 1 0 0 0]  # Assume no bias and unit rotation
    satellite_estimate = SATELLITE(_̂J, m, d, state_est, initial_covariance) # Jhat = J for now

    return satellite_truth, satellite_estimate
end

function get_random_quat()
    r = randn(3); r = r / norm(r) # random axis 
    θ = deg2rad(rand(0:360))

    q = [r * sin(θ/2); cos(θ/2)]
    return q #[0, 0, 0, 1]
end

function get_random_epoch()
    year = rand(2019:2021)
    month = rand(1:12)
    if month == 2
        day = rand(1:28)
    elseif ((month == 4) || (month == 6) || (month == 9) || (month == 11))
        day = rand(1:30)
    else
        day = rand(1:31)
    end
    hour = rand(1:23)
    min = rand(1:59)
    sec = rand(1:59)
    nano = rand(1:(1e9 - 1))

    epc_jl = Epoch(year, month, day, hour, min, sec, nano)
    epc_py = py"Epoch.ymd"(year, month, day, hour, min, sec, nano)

    return epc_jl, epc_py
end


__init_python__()

refl_dict = matread("../refl.mat") 
model = py"EarthAlbedo"(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])
refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])
cell_centers = model.get_albedo_cell_centers()
albedo = ALBEDO(refl, cell_centers)
dt = 1.0


sat, _ = generate_satellite()
obj_py = py"Simulator"(sat, dt, model)
T_py = obj_py.get_mag_calib_matrix(sat)
T_jl = generate_mag_calib_matrix(sat)

# @testset "Noise Matrix" begin
#     # No real way to compare the random things...
# end


if false 
    @testset "Mag Matrix" begin
        for i = 1:50
            sat_truth, _ = generate_satellite()

            obj_py = py"Simulator"(sat_truth, dt, model)

            T_py = obj_py.get_mag_calib_matrix(sat_truth)
            T_jl = generate_mag_calib_matrix(sat_truth)

            @test T_py ≈ T_jl 
        end 
    end

    @testset "Gyro Meas" begin
        for i = 1:50
            sat_truth, _ = generate_satellite()
            obj_py = py"Simulator"(sat_truth, dt, model)
            obj_py.gyro_noise_factor = 0.02

            x = 10 * randn(16)

            w_jl, wm_jl, η_jl = generate_gyro_measurement(x)
            w_py, wm_py, η_py =   obj_py.gyro_measurement(x)

            @test w_jl ≈ w_py
        end
    end

    @testset "GPS Meas" begin
        for i = 1:50
            sat_truth, _ = generate_satellite() 
            obj_py = py"Simulator"(sat_truth, dt, model)
            obj_py.gps_noise_factor = (5e4)

            x = 1000000 * randn(16)

            p_jl, η_jl = generate_gps_measurement(x)
            p_py, η_py = obj_py.gps_measurement(x)

            @test (p_jl-η_jl) ≈ (p_py - η_py)
        end
    end

    @testset "Diode Measurement" begin
        constants = (temp = 5, _E_am0 = 1366.9 )
        for i = 1:20
            sat, _ = generate_satellite()
            obj_py = py"Simulator"(sat, dt, model)
            obj_py.current_noise_sigma = 0.1

            pos = 1e7 * randn(3)
            t, _ = get_random_epoch()
            si  = sun_position(t)
            sb_unit = randn(3); sb_unit = sb_unit / norm(sb_unit)

            ecl = rand() > 0.5 ? 1.0 : 0

            I_py, Ĩ_py, η_py = obj_py.diode_measurement(pos, si, sb_unit, ecl)
            I_jl, Ĩ_jl, η_jl = generate_diode_currents(sat, pos, albedo, si, sb_unit, ecl, constants)

            @test I_py ≈ I_jl
        end
    end

    @testset "Sun Vectors" begin
        for i = 1:75
            sat, _ = generate_satellite()
            obj_py = py"Simulator"(sat, dt, model)
            obj_py.sun_noise_sigma = 0.0 * deg2rad(2.0)

            bRi = dcm_from_q(get_random_quat())'
            t, _ = get_random_epoch()
            pos = 1e8 * randn(3)
            if norm(pos) < R_EARTH
                continue
            end

            si_jl, sb_jl, ecl_jl = update_sun_vectors(pos, t, bRi, dt)
            si_py, sb_py, ecl_py = obj_py.update_sun_vectors(pos, t, bRi)

            @test si_jl  ≈ si_py 
            # @test sb_jl  ≈ sb_py # Requires no noise on jl and py
            @test ecl_jl ≈ ecl_py
        end
    end

end

    # @testset "Magnetometer Measurement" begin
    #     for i = 1:1
            sat, _ = generate_satellite()
            obj_py = py"Simulator"(sat, dt, model)
            obj_py.mag_noise_sigma = deg2rad(2.0)
            # pos = rand(Normal(1e7, 10000), 3)
            pos = [2.0161950101103755e6, 4.002306585404176e6, -5.111198784143326e6]
            # t, t_py = get_random_epoch()
            t = Epoch(2021, 5, 20, 12, 0, 0, 0)
            bRi = I(3) #dcm_from_q(get_random_quat())'
            
            bi_py, bb_py, b̃_py = obj_py.magnetometer_measurement(pos, t, bRi)
            bi_jl, bb_jl, b̃_jl = generate_magnetic_field(pos, t, sat, bRi, dt)

            @test bi_jl ≈ bi_py
    #     end
    # end



# ##### FULL SET #####
# @testset "Full set" begin
#     constants = (temp = 5, _E_am0 = 1366.9 )
#     sat, _ = generate_satellite()
#     obj_py = py"Simulator"(sat, dt, model)
#     obj_py.gyro_noise_factor = 0.02
#     obj_py.gps_noise_factor  = 5e4
#     obj_py.current_noise_sigma = 0.1
#     obj_py.sun_noise_sigma = 0.0 * deg2rad(2.0)

#     x = randn(16)
#     x[1:3] = 1e7 * randn(3)  # VALIDATE
#     t, t_py = get_random_epoch()

#     truth_py, sensors_py, noise_py, ecl_py = obj_py.generate_measurements(x, t_py, CONSTANTS, dt)
#     truth, sensors, ecl, noise = generate_measurements_alt(SIM(1.0), sat, albedo, x, t, CONSTANTS, dt)

#     @test ecl ≈ ecl_py 
#     @test truth.Bᴵ_hist ≈ truth_py.Bi_hist 
#     @test truth.sᴵ_hist ≈ truth_py.si_hist 
#     @test truth.Bᴮ_hist ≈ truth_py.Bb_hist 
#     @test truth.ŝᴮ_hist ≈ truth_py.sb_hist

# end


 