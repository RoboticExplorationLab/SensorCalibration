# [src/MissionSim/Simulator/Simulator.jl]

""" To Do:

    - Dynamics: see TODO in dynamics.jl 
    - Measurements: see TODO in measurements.jl
"""

module Simulator 

using LinearAlgebra, StaticArrays 
using EarthAlbedo, SatelliteDynamics 
using Distributions, Plots

include("../CustomStructs.jl"); using .CustomStructs 
include("../quaternions.jl")
include("dynamics.jl")
include("measurements.jl")
include("../mag_field.jl")


export rk4 
export generate_measurements

end