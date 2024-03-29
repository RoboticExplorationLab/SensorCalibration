# [src/Estimator/Estimator.jl]


module Estimator 

using StaticArrays, Plots, LinearAlgebra, Distributions, ForwardDiff
using SatelliteDynamics, EarthAlbedo

using Infiltrator

# Because of custom structs, I need to NOT define them again if they have been defined by a different module, but i DO need to define them if they havent been defined
if !(@isdefined STATE)
    try 
        using ..CustomStructs 
    catch 
        @info "Defining CustomStruct in Estimator..."
        include("../CustomStructs.jl"); using .CustomStructs 
    end
end

include("../quaternions.jl")
include("../mag_field.jl")


include("magnetometer_calibration.jl")
include("mekf.jl"); @info "Using FD in prediction!"

export estimate, reset_cov!
export MEKF_DATA, MAG_CALIBRATOR

end

