# [src/MissionSim/Controller/Controller.jl]

""" To Do:

    - Finish and test B-Dot controller 
    - Verify optimal κ (Reread paper...?)

    - Optimal κ; do I compute it each time? Depends on orbital period, 
        inertia, and ≈inclination, so it should vary. But I don't track 
        those things in ctrl, nor do we want to recompute every time (calc
        in config, declare as a const global/param struct...?)
"""

module Controller

# include("../CustomStructs.jl"); using .CustomStructs
using ..CustomStructs
using LinearAlgebra 
using StaticArrays

include("detumbler.jl")

export generate_command   # Interface function that is defined on all controller types
export DETUMBLER          # Controller that detumbles the CubeSat
export b_cross

end