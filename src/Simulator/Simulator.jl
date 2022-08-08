# [src/Simulator/Simulator.jl]

""" To Do:

    - Dynamics: see TODO in dynamics.jl 
    - Measurements: see TODO in measurements.jl
"""

module Simulator 

using LinearAlgebra, StaticArrays 
using EarthAlbedo, SatelliteDynamics 
using Distributions, Plots
using Infiltrator

# Because of custom structs, I need to NOT define them again if they have been defined by a different module, but i DO need to define them if they havent been defined
if !(@isdefined STATE)
    try 
        using ..CustomStructs 
    catch 
        @info "Defining CustomStruct in Simulator..."
        include("../CustomStructs.jl"); using .CustomStructs 
    end
end

include("../quaternions.jl")
include("dynamics.jl")
include("measurements.jl")
include("../mag_field.jl")


export rk4 
export generate_measurements

end