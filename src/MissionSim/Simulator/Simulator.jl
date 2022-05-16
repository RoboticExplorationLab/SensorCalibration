# [src/MissionSim/Simulator/Simulator.jl]

""" To Do:

    - Separate Dynamics, Measurements 
    - Comment, test, and speed up each 
    - Update README
    - Dynamics:
        - Add in drag that is dependent on orientation?
        - Get accurate gyro bias values
        - See if eclipse is updated (update here and in measurement)
"""

module Simulator 

include("../CustomStructs.jl"); using .CustomStructs 
include("../../../test/quaternions.jl")
using LinearAlgebra, StaticArrays 
using EarthAlbedo, SatelliteDynamics 
using Distributions

include("dynamics.jl")
include("measurements.jl")

export rk4 
export generate_measurements

end