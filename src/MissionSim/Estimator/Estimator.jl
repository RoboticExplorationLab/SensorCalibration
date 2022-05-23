# [src/MissionSim/Estimator/Estimator.jl]

""" To Do:

    - Estimator: Code, Comment, Clean, Test 
    - mag calib: Code, Comment, Clean, Test
    - diode calib: Code, Comment, Clean, Test
    - mekf: Code, Comment, Clean, Test
"""

module Estimator 

using StaticArrays, Plots, LinearAlgebra, Distributions, ForwardDiff
using SatelliteDynamics, EarthAlbedo

include("../CustomStructs.jl");  using .CustomStructs 
include("../quaternions.jl")
include("../mag_field.jl")


include("magnetometer_calibration.jl")
include("mekf.jl")

export estimate, MEKF_DATA, MAG_CALIBRATOR


end

