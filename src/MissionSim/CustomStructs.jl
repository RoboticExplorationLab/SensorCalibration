module CustomStructs
 

""" TO Do
 - Decide on default photodiode configuration (in diodes and sat state)
 - Add in plots, magnetometer noise to the 'NOISE' struct...?
 - Alphabetize? Organize somehow (in separate files?)
 - flags

 - SATELLITE has diodes in SAT_STATE && in DIODES and that can be way problematic (make 'diodes' a function call that pulls from state?)
 - SAT_STATE and DIODES both have the C, α, ϵ -> remove from sat_state


 try reduce(hcat, v) instead of (hcat([v;]...)) in plotting!

    Make a 'Structs/types' folder and have each one in its own, so we can add more functions cleanly
"""

using EarthAlbedo
using SatelliteDynamics
using StaticArrays, Distributions, LinearAlgebra
using Plots

# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo 
using Infiltrator

export ALBEDO, ESTIMATES, TRIVIAL


include("CustomStructs/diodes.jl"); export DIODES 
include("CustomStructs/flags.jl");  export FLAGS
include("CustomStructs/ground_truth.jl"); export GROUND_TRUTH 
include("CustomStructs/magnetometer.jl"); export MAGNETOMETER 
include("CustomStructs/noise.jl");  export NOISE
include("CustomStructs/sat_covariance.jl"); export SAT_COVARIANCE
include("CustomStructs/sat_state.jl"); export SAT_STATE, update_state
include("CustomStructs/satellite.jl"); export SATELLITE 
include("CustomStructs/sensors.jl");   export SENSORS 
include("CustomStructs/state.jl");     export STATE 



"""
    ALBEDO -> refl, cell_centers_ecef 

      Simple struct containing an EarthAlbedo.REFL struct (which contains 
    data, type, start_time, and stop_time), as well as the location of the
    center of each cell in Earth-Centered Earth-Fixed (ECEF) frame
"""
struct ALBEDO
    refl::REFL
    cell_centers_ecef::Array{<:Real, 3}
end



end # Module
