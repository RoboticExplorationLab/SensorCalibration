# [src/Controller/Controller.jl]

""" To Do:

    - Finish and test B-Dot controller 

    - Optimal κ; do I compute it each time? Depends on orbital period, 
        inertia, and ≈inclination, so it should vary. But I don't track 
        those things in ctrl, nor do we want to recompute every time (calc
        in config, declare as a const global/param struct...?)
"""

module Controller

# Because of the custom structs I am using, I need to NOT define them again if they have been defined by a different module, but i DO need to define them if they havent been defined
if !(@isdefined STATE)
    try 
        using ..CustomStructs 
    catch 
        @info "Defining CustomStruct in Controller..."
        include("../CustomStructs.jl"); using .CustomStructs 
    end
end


using LinearAlgebra 
using StaticArrays

include("detumbler.jl")

export generate_command   # Interface function that is defined on all controller types
export DETUMBLER          # Controller that detumbles the CubeSat
export b_cross



end