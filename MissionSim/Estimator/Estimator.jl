module Estimator

# I can just include a file for each estimator to clean up, right? 

using ..CustomStructs

using LinearAlgebra
using ForwardDiff
using SatelliteDynamics
using EarthAlbedo
# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo 
using Random, Distributions, StaticArrays, PyCall


# Primary Functions
export estimate_vals
export initialize_diodes
export initialize_magnetometer
export initialize
export save_global_variables ## REMOVE
export new_diode_calib
export new_mekf_data

export check_if_run # # # # #

# Types of estimators available:
export MAG_CALIB
export DIODE_CALIB
export MEKF_DATA

# include("estimator_config_file.jl")
const _E_am0 = 1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 
include("../rotationFunctions.jl")    # Contains general functions for working with quaternions
include("magnetometer_calibration.jl")
# include("diode_calibration.jl"); @info "Using standard MEKF for diode cal"
include("diode_calibration_sqrt.jl"); @info "Using SQRT KF"

# include("mekf.jl"); @info "Using Standard MEKF"
# include("mekf_with_consideration.jl"); @info "Using Consider MEKF, no noise!"
include("mekf_sqrt.jl"); @info "Using SQRT MEKF"

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
    diode_albedo = 0.0
    r_g = zeros(Float64, 3)
    for r = 1:1:size(albedo_matrix, 1)
        for c = 1:1:size(albedo_matrix, 2)
            if albedo_matrix[r,c] != 0
                r_g .= view(cell_centers_ecef, r, c, :) .- sat_pos
                r_g .= r_g ./ norm(r_g)  # Make unit

                cell_albedo = (albedo_matrix[r,c] * dot(surface_normal, r_g))

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo += cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end

function initialize_magnetometer()
    return MAGNETOMETER(ones(3), zeros(3), zeros(3)) # Assume no bias, no non-orthogonality, and a unit scale factor 
end

function initialize_diodes(_num_diodes)
    """ Initial estimate for photodiode parameters (what they were intended to be when installed) """
    
    scale_factor_init_est = ones(_num_diodes) 
    azi_angles_init_est   = [0.0;      pi;    (pi/2); (-pi/2);  0.0;    pi] 
    elev_angles_init_est  = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)]
    
    diode = DIODES(scale_factor_init_est, azi_angles_init_est, elev_angles_init_est)

    return diode
end

####################################################################
#               TRIVIAL CASE (does nothing)                        #
####################################################################

function estimate_vals(sat::SATELLITE, data::TRIVIAL)
    return sat, data 
end

end