module Estimator

using ..CustomStructs

using LinearAlgebra, ForwardDiff, SatelliteDynamics
using Random, Distributions, StaticArrays
using JLD2, MAT

using PyCall, Infiltrator, BenchmarkTools, Test

using EarthAlbedo
# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo # For testing purposes 


# Primary Functions
export estimate_vals
export initialize_diodes
export initialize_magnetometer
export new_diode_calib
export new_mekf_data

export check_if_run # # # # #

# Types of estimators available:
export MAG_CALIB
export DIODE_CALIB
export MEKF_DATA

include("estimator_config_file.jl")   # Assorted, relevant parameters
include("../rotationFunctions.jl")    # Contains general functions for working with quaternions

include("magnetometer_calibration.jl")
include("diode_calibration_sqrt.jl")
include("mekf_sqrt.jl")


####################
# SHARED FUNCTIONS #
####################

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
    """ Initializes a MAGNETOMETER object (assumes no bias, perfectly orthogonal axes, and unit scale factors """
    return MAGNETOMETER(ones(3), zeros(3), zeros(3)) 
end

function initialize_diodes(_num_diodes)
    """ Initial estimate for photodiode parameters (what they were intended to be when installed) """
    
    scale_factor_init_est = ones(_num_diodes) 
    azi_angles_init_est   = [0.0;      pi;    (pi/2); (-pi/2);  0.0;    pi] 
    elev_angles_init_est  = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)]
    
    diode = DIODES(scale_factor_init_est, azi_angles_init_est, elev_angles_init_est)

    return diode
end

function generate_mag_calib_matrix(sat::SATELLITE)
    """ Generates the calibration matrix that alters the measured magnetic field vector in body frame """
    a, b, c = sat.magnetometer.scale_factors
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles

    T = [a        0.0              0.0;
        b*sin(ρ)  b*cos(ρ)         0.0;
        c*sin(λ)  c*sin(ϕ)*cos(λ)  c*cos(ϕ)*cos(λ)]

    return T
end


####################################################################
#               TRIVIAL CASE (does nothing)                        #
####################################################################

function estimate_vals(sat::SATELLITE, data::TRIVIAL)
    return sat, data 
end

end