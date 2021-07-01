module Estimator

# I can just include a file for each estimator to clean up, right? 

using ..CustomStructs

using LinearAlgebra
using ForwardDiff
using SatelliteDynamics, EarthAlbedo  # Used for sun position in MEKF
using Random, Distributions


# Primary Functions
export estimate_vals
export initialize

# Types of estimators available:
export MAG_CALIB
export DIODE_CALIB
export MEKF


include("estimator_config_file.jl")
include("../rotationFunctions.jl")    # Contains general functions for working with quaternions
include("magnetometer_calibration.jl")
include("diode_calibration.jl")
include("mekf.jl")


function compute_diode_albedo(albedo_matrix, cell_centers_ecef, surface_normal, sat_pos)
    """ 
        Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
            = cell_albedo * surface_normal^T * r_g, with r_g as a unit vector in the direction of the grid point on Earth

        Arguments:
        - albedo_matrix: Albedo values for each cell on the Earth's surface         | [num_lat x num_lon] 
        - surface_normal: Photodiode surface normal                                 | [3,]
        - sat_pos: Cartesian position of satellite                                  | [3,]

        Returns:
        - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
    """    
    cell_albedos = zeros(size(albedo_matrix))

    rows, cols = size(albedo_matrix)
    diode_albedo = 0.0
    for r = 1:1:rows
        for c = 1:1:cols
            if albedo_matrix[r,c] != 0
                r_g = cell_centers_ecef[r,c,:] - sat_pos # Distance from satellite to cell center
                r_g = r_g / norm(r_g)  # Make unit

                cell_albedo = (albedo_matrix[r,c] * (surface_normal * r_g))[1]

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo = diode_albedo + cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end

function initialize(m::MAGNETOMETER)
    return m 
end

# MOVE to CONFIG? σᶜ, μᶜ, etc...?
function initialize(sat::SATELLITE, diode::DIODES)
    """ SAT TURTH not EST"""
    _num_diodes = length(sat.diodes.azi_angles)
    
    # # Calib Azi Elev ∈ DIODE 
    # scale_factor_init_est = rand(Normal(1.0, 0.1), _num_diodes)
    # azi_angles_init_est   = rand(Normal(pi, deg2rad(7)), _num_diodes)       
    # elev_angles_init_est  = rand(Normal(pi/2, deg2rad(7)), _num_diodes)
    # diode = DIODES(scale_factor_init_est, azi_angles_init_est, elev_angles_init_est)


    scale_factor_init_est = sat.diodes.calib_values + rand(Normal(0.0, 0.3), _num_diodes)
    azi_angles_init_est   = sat.diodes.azi_angles + rand(Normal(0.0, deg2rad(10)), _num_diodes)       
    elev_angles_init_est  = sat.diodes.elev_angles + rand(Normal(0.0, deg2rad(10)), _num_diodes)
    diode = DIODES(scale_factor_init_est, azi_angles_init_est, elev_angles_init_est)


    return diode
end

####################################################################
#               TRIVIAL CASE (does nothing)                        #
####################################################################

function estimate_vals(sat::SATELLITE, data::TRIVIAL)
    return sat, data, false 
end





end