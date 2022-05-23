# [src/MissionSim/Estimator/Estimator.jl]

""" To Do:

    - Estimator: Code, Comment, Clean, Test 
    - mag calib: Code, Comment, Clean, Test
    - diode calib: Code, Comment, Clean, Test
    - mekf: Code, Comment, Clean, Test
"""

module Estimator 

using StaticArrays, Plots, LinearAlgebra, Distributions, ForwardDiff

include("../CustomStructs.jl");  using .CustomStructs 


include("magnetometer_calibration.jl")
include("diode_calibration.jl")
include("mekf.jl")

export estimate





end

