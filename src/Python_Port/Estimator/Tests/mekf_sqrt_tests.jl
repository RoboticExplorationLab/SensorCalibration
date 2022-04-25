"""
    TEST FILE for comparing mekf_sqrt.jl with sqrtSatMekf.py 
        to verify the python implementation.

    Data is specified and then specific functions common to both files are run.
    A full iteration of the MEKF is also run and compared
"""

using Infiltrator, PyCall, Test, MAT 
using SatelliteDynamics, LinearAlgebra, EarthAlbedo 

include("../mag_field.jl")            # Contains IGRF13 stuff
include("../rotationFunctions.jl");   # Contains general functions for working with attitude and quaternions
include("../CustomStructs.jl");         using .CustomStructs 
include("../system_config_file.jl"); 
include("Estimator.jl");   using .Estimator

function __init_python__()
    # Allows us to run python scripts with PyCall
    py"""
    import sys 
    sys.path.insert(0, "./") # Allows us to include nearby python scripts

    import numpy as np 

    from earthAlbedo import EarthAlbedo 
    from satelliteDynamics import *
    from sqrtSatMekf import SqrtSatMekf

    import scipy.io
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


### Set up test data 
__init_python__()

refl_dict = matread("../refl.mat") 
model = py"EarthAlbedo"(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])
refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])
cell_centers = model.get_albedo_cell_centers()
albedo = ALBEDO(refl, cell_centers)

if true
    @testset "Individual Functions" begin 

        ### Set up necessary data
        sensors = SENSORS( [20.0, 7.0, 10.0], [1.0, 0.0, 0.67, 0.33, 0.1, 0.9],
                        [1.2, 4.4, 1.1], [10000, 10000, 10000000])

        q = [0, 0, 0, 1]

        sat_true, sat_est = generate_satellite()

        num_diodes = 6
        σ_q = (10*pi/180) 
        σ_β = (10*pi/180)
        p = [σ_q * ones(3); σ_β * ones(3)].^2
        P₀ = diagm(p)
        P₀ = chol(P₀)
        sat_est.covariance = P₀


        est = new_diode_calib(albedo, sensors, SYSTEM, q, sat_est)
        sᴵ  = [0.0, 0.0, 1e10] - sensors.gps
        Bᴵ  = [33, 2.0, 11.3]
        sᴮ  = [0.0, 0.0, 1.0]
        Bᴮ  = sensors.magnetometer

        est.inertial_vecs = [sᴵ Bᴵ]'
        est.body_vecs = [sᴮ Bᴮ]'
        est.current_meas = sensors.diodes

        t = Epoch(2021, 9, 1, 11, 0, 0, 0.0)
        t_py = py"Epoch.ymd"(2021, 9, 1, 11, 0, 0, 0.0)

        x = est.sat_state
        Mm = randn(5,5); Mm = Mm * Mm'

        # UPDATE JULIA 
        jl_data = new_diode_calib(albedo, sensors, SYSTEM, q, sat_true)
        jl_data.sat_state = est.sat_state 
        jl_data.ang_vel = sensors.gyro 
        jl_data.pos = sensors.gps
        jl_data.current_meas = sensors.diodes 
        jl_data.inertial_vecs = est.inertial_vecs
        jl_data.body_vecs = est.body_vecs
        jl_data.time = t 
        
        # UPDATE PYTHON 
        data = Dict("time" => t_py, "q" => est.sat_state[1:4], "bias" => est.sat_state[5:7],
                    "sun_inert" => est.inertial_vecs[1,:], "mag_inert" => est.inertial_vecs[2,:],
                    "sun_body" => est.body_vecs[1,:], "mag_body" => est.body_vecs[2,:],
                    "diodes" => sensors.diodes, "gyro" => sensors.gyro, "gps" => sensors.gps)
        
        c = sat_true.diodes.calib_values    
        a = sat_true.diodes.azi_angles 
        e = sat_true.diodes.elev_angles 

        ### Begin Tests

        # NOISE MATRICES 
        est = new_diode_calib(albedo, sensors, SYSTEM, q, sat_true)
        obj_py = py"SqrtSatMekf"(c, a, e, SYSTEM, model, data)
        @test est.W[1:6, 1:6] ≈ obj_py.W 
        @test est.V ≈ obj_py.V


        # CHOLESKY
        obj_py = py"SqrtSatMekf"(c, a, e, SYSTEM, model, data)
        chol_py = obj_py.chol(Mm)
        chol_jl = chol(Mm)
        @test chol_py ≈ chol_jl

        # PREDICTION
        obj_py = py"SqrtSatMekf"(c, a, e, SYSTEM, model, data)
        xp_py, Pp_py = obj_py.predict()
        xp_jl, A_jl = prediction(x, sensors.gyro, SYSTEM._dt)
        @test xp_py ≈ xp_jl

        # MEASUREMENT 
            # MAGNETOMETER MEASUREMENT 
        obj_py = py"SqrtSatMekf"(c, a, e, SYSTEM, model, data)
        Bi_unit = data["mag_inert"] / norm(data["mag_inert"])
        ymag_py, Cmag_py = obj_py.mag_measurement(est.sat_state, Bi_unit)
        ymag_jl, Cmag_jl = mag_measurement(est.sat_state, Bi_unit)
        @test ymag_py ≈ ymag_jl 
        @test Cmag_py ≈ Cmag_jl 

            # CURRENT MEASUREMENT 
        obj_py = py"SqrtSatMekf"(c, a, e, SYSTEM, model, data)
        si_unit = data["sun_inert"] / norm(data["sun_inert"])
        ymag_py, Cmag_py = obj_py.current_measurement(est.sat_state, si_unit)
        ymag_jl, Cmag_jl = current_measurement(est.sat_state, c, a, e, 
                                                si_unit, SYSTEM._num_diodes, sensors.gps,
                                                t, albedo )
        @test ymag_py ≈ ymag_jl
        @test Cmag_py ≈ Cmag_jl 
    end
end



### HELPER FUNCTIONS FOR TESTING 
function get_random_quat()
    r = randn(3); r = r / norm(r) # random axis 
    θ = deg2rad(rand(0:360))

    q = [r * sin(θ/2); cos(θ/2)]
    return q #[0, 0, 0, 1]
end

function fill_sensors()
    mag_meas = abs.(10 * randn(3)) #[20.0, 7.0, 10.0]

    d = 0.5 .+ 0.25 * randn(3)    
    diode_meas = abs.([ d[1], 1.0-d[1], d[2], 1.0 - d[2], d[3], 1 - d[3] ])  # [1.0, 0.0, 0.67, 0.33, 0.1, 0.9]

    gyro_meas = 1 .+ randn(3) # [1.2, 4.4, 1.1]

    gps_meas = 0
    while norm(gps_meas) < 7000e3
        gps_meas = 8000e3 * randn(3)  # [10000, 10000, 10000000]
    end

    return SENSORS(mag_meas, diode_meas, gyro_meas, gps_meas)
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

function estimate_sun_vector(sens::SENSORS, sat_est::SATELLITE)
    """ Estimates a (unit) sun vector using the diode measurements 
            (Note that the equation is strange because a 45° rotation @ ŷ was used to avoid boundary problems with the elevation angles) """

    if true # norm(sens.diodes) > eclipse_threshold  # If not eclipsed
        x₁ = (sens.diodes[1]/sat_est.diodes.calib_values[1])
        x₂ = (sens.diodes[2]/sat_est.diodes.calib_values[2])
        y₁ = (sens.diodes[3]/sat_est.diodes.calib_values[3])
        y₂ = (sens.diodes[4]/sat_est.diodes.calib_values[4])
        z₁ = (sens.diodes[5]/sat_est.diodes.calib_values[5])
        z₂ = (sens.diodes[6]/sat_est.diodes.calib_values[6])

        sun_vec_est = [(x₁*cos(-pi/4) + z₁*cos(pi/4) + x₂*cos(3*pi/4) + z₂*cos(-3*pi/4));    
                        y₁ - y₂;
                       (x₁*cos(3*pi/4) + z₁*cos(pi/4) + x₂*cos(-pi/4) + z₂*cos(-3*pi/4))] 

        sun_vec_est /= norm(sun_vec_est)
    else
        # Assume that we are in eclipse
        sun_vec_est = [0; 0; 0]
    end
        
    return sun_vec_est # Unit - ŝᴮ
end

@testset "Multiple Iterations" begin 
    for i = 1:10
        sat_true, sat_est = generate_satellite()

        num_diodes = 6 
        σ_q = (10*pi/180) 
        σ_β = (10*pi/180)

        p = [σ_q * ones(3); σ_β * ones(3)].^2
        P₀ = diagm(p)
        P₀ = chol(P₀)
        sat_est.covariance = P₀

        sensors = fill_sensors()
        q = get_random_quat() 
        R = dcm_from_q(q)' 

        sᴮ = estimate_sun_vector(sensors, sat_est)
        sᴵ = 1e10 * (R' * sᴮ) - sensors.gps
        Bᴮ  = sensors.magnetometer
        Bᴵ  = R' * Bᴮ

        est = new_diode_calib(albedo, sensors, SYSTEM, q, sat_est)
        est.inertial_vecs = [sᴵ Bᴵ]'
        est.body_vecs = [sᴮ Bᴮ]'
        est.current_meas = sensors.diodes

        t, t_py = get_random_epoch()

        ### JULIA (REFERENCE) IMPLEMENTATION 
        jl_data = new_mekf_data(albedo, sensors, SYSTEM, q, sat_true)
        jl_data.sat_state = est.sat_state 
        jl_data.ang_vel = sensors.gyro 
        jl_data.pos = sensors.gps
        jl_data.current_meas = sensors.diodes 
        jl_data.inertial_vecs = est.inertial_vecs
        jl_data.body_vecs = est.body_vecs
        jl_data.time = t 
        jl_sat, data = estimate_vals(sat_true, jl_data)

        ### PYTHON (NEW) IMPLEMENTATION 
        # PYTHON METHOD
        data = Dict("time" => t_py, "q" => est.sat_state[1:4], "bias" => est.sat_state[5:7],
                    "sun_inert" => est.inertial_vecs[1,:], "mag_inert" => est.inertial_vecs[2,:],
                    "sun_body" => est.body_vecs[1,:], "mag_body" => est.body_vecs[2,:],
                    "diodes" => sensors.diodes, "gyro" => sensors.gyro, "gps" => sensors.gps)

  
        c = sat_true.diodes.calib_values    
        a = sat_true.diodes.azi_angles 
        e = sat_true.diodes.elev_angles 
                
        obj_py = py"SqrtSatMekf"(c, a, e, SYSTEM, model, data)
        _ = obj_py.update_estimate(data, sat_est)

        
        ### COMPARISON 
        @test obj_py.q ≈ jl_sat.state[1:4]
        @test obj_py.bias ≈ jl_sat.state[5:7]
    end
end





