module Simulator
# Interface with Dynamics and Measurements files

using LinearAlgebra, SatelliteDynamics 
using Random, Distributions 
using Plots, SparseArrays

using EarthAlbedo
# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo;    # ← Used for testing new code

using BenchmarkTools, Infiltrator
# using PyCall

using ..CustomStructs
include("../rotationFunctions.jl"); # Helpful functions for working with quaternions
include("../mag_field.jl")  # Used for mag field values

include("simulator_config_file.jl")
include("Dynamics.jl")
include("Measurements.jl")

export run_dynamics 
export generate_measurements


const with_noise = true    # Used for debugging
if (!with_noise)
    @info "No noise in Simulator!"
end


# TODO can I just delete this...
export SIM
struct SIM # For now, just used for multiple dispatch 
    junk
end

end