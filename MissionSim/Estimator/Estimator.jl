module Estimator

# I can just include a file for each estimator to clean up, right? 

include("estimator_config_file.jl")
using ..CustomStructs

using LinearAlgebra
using ForwardDiff
using SatelliteDynamics, EarthAlbedo  # Used for sun position in MEKF
using Random, Distributions


export estimate_vals
export initialize

# Types of estimators available:
export MAG_CALIB
export DIODE_CALIB
export MEKF

include("../rotationFunctions.jl")    # Contains general functions for working with quaternions

include("magnetometer_calibration.jl")
include("diode_calibration.jl")




function initialize(m::MAGNETOMETER)
    return m 
end


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

# MAKE SUBMODULES...?

####################################################################
#               TRIVIAL CASE (does nothing)                        #
####################################################################
# struct TRIVIAL
#     junk
# end

function estimate_vals(sat::SATELLITE, data::TRIVIAL)
    t = data
    return sat, false 
end







####################################################################
#                      STANDARD MEKF                               #
####################################################################
struct MEKF 
    # state 
    # covariance 
    # inertial_vecs 
    # ang_vel 
    # body_vecs 
    # W 
    # V

    dt 
end

function estimate_vals(sat::SATELLITE, data::MEKF)
    # Initialize with end of diode_calib
    return sat, true
end

function initialize(data::MEKF)
    # In place
    return data
end


end